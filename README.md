# Simulating-New-Keynesian-Model-with-Cost-Push-Shock-using-JuliaPerturbation-Toolbox
This work focuses on a medium-scale New-Keynesian model in a closed economy with com-
plete financial markets and capital accumulation, and its implementation in Julia using the
JuliaPerturbation toolbox created by Alvaro Salazar-Perez and Hernan D. Seoane (2023).
The model is extended by adding a cost-push shock, which is described in detail in a paper
by Peter N. Ireland (2002) . The impulse response functions (IRFs) are simulated for three
different types of shocks: a monetary policy shock, a TFP shock, and a cost-push shock.
The JuliaPerturbation toolbox is used to perform the perturbation and estimation of the
model, and finally the Impulse Response Functions are compared and discussed.

In the zip folder you can find the following documents:

 1. solution_functions_7_0.jl: This file is from Salazar and Seoane Julia Toolbox, that is available in: https://github.com/HernanSeo/JuliaPerturbation
 
 2. Impulse Response Analysis.jl: In this file I solve the model, and plot the IRF for three different identification shocks 

 3. Determinancy Region Analysis.jl: In this file I study indeterminancy for α_y and α_ϕ for fixed intervals, iteratively solving the model and storing the values of the parameters that lead to a stable solution.

 4. Figures: A folder with the output of the previous code

 5. Writing Output: The final document in which the results and simulations are analysed in a research paper fashion. 

 ### NOTE!: To properly understand the exercise please read the Writing Output file ### 
