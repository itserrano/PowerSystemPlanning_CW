# Power System Planning Coursework
Github repository for the Power System Planning Final Project

### Instructions
- Folder __Data__ contains the used data for the model. All the values used are obtained from the references presented in the References section.
- The model is present in the file `PSP_CW.jl`. This file requires the following Julia packages to run:
    > `JuMP`, `DataFrames`, `PlotlyJS`, `CSV`, `XLSX`, `Statistics`, `DataStructures`, `Dates`, `Gurobi`
- Furthermore, additional files for auxiliary functions are provided. `Data_Functions.jl` have functions to read the data and prepare it for the optimisation model, and `mycolorscheme.jl` defines a new colour palette for the plots. These files require the following packages:
    > `Colors`, `ColorSchemes`
- Results of the problem are saved in the folder __Results__. The output is one Excel file per scenario, with a generation sheet, installed capacity sheet and emissions sheet.


### References
- *NREL (National Renewable Energy Laboratory)*. 2023. "2023 Annual Technology Baseline." Golden, CO: National Renewable Energy Laboratory. [https://atb.nrel.gov/](https://atb.nrel.gov/). 
- *Renewables Ninja*. 2023. Wind and Solar Data for United Kingdom and Spain. [https://www.renewables.ninja/](https://www.renewables.ninja/)
- Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2023) – with major processing by Our World in Data. “Coal” [dataset]. Energy Institute, “Statistical Review of World Energy” [original data].
- Energy Institute based on S&P Global Platts - Statistical Review of World Energy (2023) – with major processing by Our World in Data. “Gas price” [dataset]. Energy Institute, “Statistical Review of World Energy” [original data].
- National Grid ESO. 2024. "FES 2023 Data workbook". [https://www.nationalgrideso.com/document/283061/download](https://www.nationalgrideso.com/document/283061/download).
- European Commission. 2023. "EU energy statistical pocketbook and country datasheets". [https://energy.ec.europa.eu/data-and-analysis/eu-energy-statistical-pocketbook-and-country-datasheets_en](https://energy.ec.europa.eu/data-and-analysis/eu-energy-statistical-pocketbook-and-country-datasheets_en)