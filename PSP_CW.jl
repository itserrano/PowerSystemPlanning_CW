using JuMP,DataFrames, CSV, Gurobi, Statistics
import ParametricOptInterface as POI
include("Data_Functions.jl")
root = "Data/"

# Load data
gen_data, gen_set = generator_data("Data/CostsData.csv");
demand_df = CSV.read("Data/demanddata_2022.csv", DataFrame);
D, date = demand_df.ND, demand_df.Date;
H = 1:length(D);
VoLL = 9000;
cf = capacity_factors("Data/TMY.csv");
current_set = current_capacity("Data/UK_GenMix_2022.csv");
tx_OrderedDict = transmission_costs("Data/LineCosts.csv");
import_price = CSV.read("Data/EU-Spain.csv",DataFrame)[!,"Price (EUR/MWhe)"];

G, G_rve, G_bess, G_EU, G_conv = gen_set
existing, cap_0, tau = current_set

carbon_tax, CBA, emissions = 0.0 , 0.0 ,zeros(8760,1)

# Define model
function Expansion_Model()
    global carbon_tax, CBA, emissions
    global cap, gen, nse, soe, charge, discharge, rep_gen, flow, line_cap

    ExpansionModel = Model(Gurobi.Optimizer)
    set_optimizer_attribute(ExpansionModel, "Threads", 16)
    @variables(ExpansionModel, begin
        cap[g in G] >=0
        gen[g in [G_conv;G_rve;G_EU], t in H] >=0
        nse[t in H]>=0
        soe[g in G_bess, t in H] >= 0
        charge[g in G_bess, t in H] >= 0
        discharge[g in G_bess, t in H] >= 0
        rep_gen[g in G_EU, t in H] >= 0
        flow[t in H] >= 0
        line_cap >= 0
    end)

    @constraints(ExpansionModel, begin
        EU_flow[t in H], flow[t] == sum(gen[g,t] for g in G_EU)
        EU_gen[t in H], sum(gen[g,t] for g in G_EU) == sum(rep_gen[g,t] for g in G_EU)
        demandSupply[t in H], sum(gen[g,t] for g in [G_conv;G_rve]) + sum((discharge[g,t]-charge[g,t]) for g in G_bess) + nse[t] + flow[t] == 1.6*D[t]

        lineConstraint[t in H], flow[t] <= line_cap + tx_OrderedDict["CAPACITY"]

        capacityConstraintRen[g in G_rve, t in H], gen[g,t] <= (cap[g]+cap_0[g]) * cf[g][t]
        capacityConstraintConv[g in G_conv, t in H], gen[g,t] <= (cap[g]+cap_0[g]) * gen_data["CapacityFactor"][g]

        capacityConstraintEU[g in G_EU, t in H], gen[g,t] <= (cap[g]+cap_0[g])*cf[g][t]
        nonFlexNuclear[t in H], gen["Nuclear",t] == (cap["Nuclear"]+cap_0["Nuclear"])*gen_data["CapacityFactor"]["Nuclear"]
        
        soeConstraint[g in G_bess, t in H; t<8760], soe[g,t+1] == soe[g,t]+(charge[g,t]*0.85-discharge[g,t]/0.85)
        soeInitial[g in G_bess, t in [1,8760]], soe[g,t] == 0.5*(cap[g]+cap_0[g])*tau[g]
        soeLimit[g in G_bess, t in H], soe[g,t] <= (cap[g]+cap_0[g])*tau[g]
        chargeLimit[g in G_bess, t in H], charge[g,t] <= (cap[g]+cap_0[g])
        dischargeLimit[g in G_bess, t in H], discharge[g,t] <= (cap[g]+cap_0[g])

        initial_discharge[g in G_bess], discharge[g,1] <= soe[g,1]
        initial_charge[g in G_bess], charge[g,1] <= (cap[g]+cap_0[g])-soe[g,1]
        
        rampUp[g in [G_conv;G_rve], t in H; 1<t], gen[g,t]-gen[g,t-1] <= gen_data["Ramping"][g]*cap[g]
        rampDown[g in [G_conv;G_rve], t in H; 1<t], gen[g,t-1]-gen[g,t] <= gen_data["Ramping"][g]*cap[g]
    end)

    @objective(ExpansionModel, Min, sum(
        (gen_data["OM"][g]+gen_data["CAPEX"][g])*(cap[g]) for g in G)+
        sum((gen_data["VarOM"][g]+gen_data["HeatRate"][g]*gen_data["FuelCost"][g])*gen[g,t] for g in [G_conv;G_rve], t in H)+
        sum(VoLL*nse[t]+flow[t]*import_price[t] for t in H) + sum(carbon_tax*gen_data["HeatRate"][g]*gen_data["CO2"][g]*gen[g,t] for g in G_conv, t in H)
        + line_cap*tx_OrderedDict["CAPEX"] + sum(carbon_tax*CBA*gen_data["HeatRate"][g]*gen_data["CO2"][g]*rep_gen[g,t] for t in H for g in ["Gas-CC-EU"])+
        sum(carbon_tax*flow[t]*emissions[t] for t in H)
    )
    optimize!(ExpansionModel)
end

@def compile_results begin 
    global useful_g = OrderedDict("Conv"=>[],"RVE"=>[],"BESS"=>[],"EU"=>[])
    for g in [G_conv;G_rve;G_EU]
        if sum(value.(gen[g,:])) > 10
            if g in G_conv
                push!(useful_g["Conv"],g)
            elseif g in G_rve
                push!(useful_g["RVE"],g)
            elseif g in G_EU
                push!(useful_g["EU"],g)
            else
                push!(useful_g["BESS"],g)
            end
        end
        println(g,": ",value(cap[g]), " MW")
    end
    for g in G_bess
        if sum(value.(discharge[g,:]) ) > 0 || sum(value.(charge[g,:])) > 0
            push!(useful_g["BESS"],g)
        end
        println(g,": ",value(cap[g]), " MW")
    end
end

scenarios = OrderedDict("Base"=>(0,0,zeros(8760,1)),
    "Carbon_Tax_UK_low"=>(0.05,0,zeros(8760,1)),
    "Carbon_Tax_UK_high"=>(0.2,0,zeros(8760,1)),
    "Carbon_Tax"=>(0.1,0,zeros(8760,1)),
    "CBA"=>(0.1,1,zeros(8760,1))
)
results = OrderedDict()
for scenario in scenarios
    name, params = scenario
    global carbon_tax, CBA, emissions = params
    Expansion_Model()
    results[name] = OrderedDict("Capacity" => value.(cap), "Generation" => value.(gen), "Flow" => value.(flow), "Line_Cap" => value(line_cap),
    "SOE" => value.(soe), "Charge" => value.(charge), "Discharge" => value.(discharge), "Rep_Gen" => value.(rep_gen), "NonServedEnergy" => value.(nse))
    @compile_results
end

# We process the data now
using DataStructures

processed_results = OrderedDict()
k = DataFrame()
for generator in G
    global k[!,generator] = [cap_0[generator]]
end

processed_results["Existing"] = OrderedDict("Installed_Capacity" => k)

for scenario in results
    name, data = scenario
    global installed_capacity = DataFrame()
    generation = DataFrame()
    for generator in G
        installed_capacity[!,generator] = [data["Capacity"][generator] + cap_0[generator]]
        if generator ∉ (G_bess)
            generation[!,generator] = data["Generation"][generator,:].data
        else
            generation[!,generator] = data["Discharge"][generator,:].data - data["Charge"][generator,:].data
        end
    end
    generation[!,"NSE"] = data["NonServedEnergy"].data
    generation[!,"Flow"] = data["Flow"].data
    for generator in G_EU
        generation[!,"Reported_"*generator] = data["Rep_Gen"][generator,:]
    end
    processed_results[name] = OrderedDict("Installed_Capacity" => installed_capacity, "Generation" => generation,"Line_Cap"=>data["Line_Cap"])
end


# Save the results in an excel file for each scenario, using CLSX package
using XLSX
for scenario in scenarios
    name, value = scenario
    value = processed_results[name]
    XLSX.writetable("Results/$(name).xlsx", "Installed_Capacity"=>value["Installed_Capacity"], "Generation"=>value["Generation"], overwrite=true)
    # XLSX.writetable("Results/$(name).xlsx", "Generation", data["Generation"])
    # XLSX.writetable("Results/$(name).xlsx", "NonServedEnergy", DataFrame(data["NonServedEnergy"]))
    # XLSX.writetable("Results/$(name).xlsx", "Flow", DataFrame(data["Flow"]))
    # XLSX.writetable("Results/$(name).xlsx", "Line_Cap", DataFrame(data["Line_Cap"]))
end

using PlotlyJS
plot(
    [bar(name=g, x=keys(processed_results), y = [processed_results[scenario]["Installed_Capacity"][1,g] for scenario in keys(processed_results)]) for g in G]
,Layout(barmode="stack"))

plot([
    bar(name="Non Served Energy", x=keys(scenarios), y = [sum(processed_results[scenario]["Generation"][:,"NSE"]) for scenario in keys(scenarios)])])

p = make_subplots(rows=length(scenarios), cols=1)
for (i,scenario) in enumerate(keys(scenarios))
    add_trace!(p, pie(labels=G, values=[sum(processed_results[scenario]["Generation"][:,g]) for g in G]), row=i, col=1)
end


for scenario in keys(scenarios)
    display(plot(
        [pie(labels=G, values=[sum(processed_results[scenario]["Generation"][:,g]) for g in G], title=scenario, textinfo="percent+label", texttemplate="%{label}    %{percent}")],
        Layout(font=attr(size=12))))
end



