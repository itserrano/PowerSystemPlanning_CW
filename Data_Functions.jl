using CSV, DataFrames, JuMP

function generator_data(file)
    generators = CSV.read(file, DataFrame)
    G = generators.Gen;
    to_use = names(generators)[.!occursin.("KW",names(generators)).&.!occursin.("Gen",names(generators))]

    G_rve = G[occursin.("Wind",G) .| (occursin.("Solar",G).& .!occursin.("EU",G))]
    G_bess = G[occursin.("BESS",G)]
    G_EU = G[occursin.("EU",G)]
    G_conv = G[G.∉ Ref([G_rve;G_bess;G_EU])]
    
    gen_data = Dict()
    for column in to_use
        colname = split(column)[1]
        gen_data[colname] = JuMP.Containers.DenseAxisArray(generators[!,column],G)
    end

    data = Nothing
    return gen_data, [G, G_rve, G_bess, G_EU, G_conv];
end

function capacity_factors(file)
    cap_names = ["OnshoreWind","OffshoreWind","Solar","Solar-EU"]
    cf_df = CSV.read(file,DataFrame)
    rename!(cf_df,["onshore","offshore","solar","solar-ES"].=>["OnshoreWind","OffshoreWind","Solar","Solar-EU"])

    cf = Dict()
    for col in cap_names
        cf[col] = JuMP.Containers.DenseAxisArray(cf_df[!,col],H)
    end
    cf["Gas-CC-EU"] = ones(8760,1);
    return cf;
end

function current_capacity(file)
    existing = CSV.read(file,DataFrame);
    cap_0 = JuMP.Containers.DenseAxisArray(existing.CAP,existing.Gen);
    G_bess = existing.Gen[occursin.("BESS",existing.Gen)]
    tau = JuMP.Containers.DenseAxisArray([2,6,10],G_bess);
    return [existing, cap_0, tau];
end

function transmission_costs(file)
    tx_costs_df = CSV.read(file,DataFrame)
    tx_dict = Dict()
    for col in names(tx_costs_df)
        tx_dict[col] = tx_costs_df[1,col]
    end
    return tx_dict;
end

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
  end

