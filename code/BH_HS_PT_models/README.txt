In this folder there is code to run the the folowing analysis

1. Simulate 1,000 biomass data sets, each over 50 time steps, for each model and several levels of environmental noise
2. Fit stock recruitment relationships to each data set using the true model and the two wrong models (the ones not chosen to generate the data)
3. Generate the optimal escapement rule given the model choice and fitted parameter values
4. Simulate that escapement rule in 1000 timesteps of fishing, where the simulation is generated using the analgous true model from 1
5. Compare the yield from each fitted model to the true optimal escapement rule

The files

'Function_Analysis.R' contains the main function that conducts the above analysis for a given set of parameters, and a list of models
'Functions_Escapement.R' contains subroutines for the 
  (a) escapement functions, 
  (b) optimal escapement formulas, 
  (c) generating the data given a parameter set
  (d) computing sum of squared error on either the log or raw data
  (e) fitting the parameters using optim with the functions in (a) and (d)
  (f) simulating an escapement strategy 
      (note routine c calls this function with effectively infinite escapement, 
      if biomass data is generated through simulation rather than randomly)
      
There are also PDFs of some results and latex files with a 0th draft of the paper
