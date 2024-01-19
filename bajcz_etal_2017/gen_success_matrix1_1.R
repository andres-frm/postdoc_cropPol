#The gen_Success_Matrix function creates a matrix that has a row and a column for every plant in the field. In the interior,
#it stores the maximum probability of successful pollination for any pollination event, such that the rows represent donors
#and columns represent receivers. Both are ordered identically such that the diagonal represents self-pollination events.

gen_Success_Matrix = function (mean_Crossing_Prob, sd_Crossing_Prob, grid_Dim, selfing_Prob_Adjust) {
  
  #The function starts by drawing a average pollen donation and reception probability for each plant.
  donate = rtruncnorm(grid_Dim^2, 0, 1, mean_Crossing_Prob, sd_Crossing_Prob)
  receive = rtruncnorm(grid_Dim^2, 0, 1, mean_Crossing_Prob, sd_Crossing_Prob)
  
  gen_Load = matrix((donate+receive)/2, grid_Dim, grid_Dim) #Store the gen_Load of each clone for reference later.
  
  #The combine function will create an average of these two probabilities for each plant and arrange them in a matrix.
  combine = function(vec1, vec2) {
    return(apply(cbind(vec1, vec2),1,mean))
  }
  combine = cmpfun(combine)
  Pol = outer(donate, receive, combine)
  
  #The final average crossing success rate (as it applies to both being a donor and a receiver) is redrawn using the average
  #calculated above as the mean and the original standard deviation. 
  for (cell in 1:length(Pol)) { 
    Pol[cell] = rtruncnorm(1, 0, 1, Pol[cell], sd_Crossing_Prob)
  }
  
  #Then, the function calculates a selfing success rate using each plant's ability to donate pollen. These selfing rates 
  #are put into the Pol matrix created earlier along the diagonal.
  selfrates = rep(0, grid_Dim^2) 
  rate.val = rep(0, grid_Dim^2) 
  for (i in 1:(grid_Dim^2)) {
    rate.val[i] = (1/((donate[i]*donate[i])+selfing_Prob_Adjust))
    if (rate.val[i] < 0) {
      rate.val[i] = 0
    }
    selfrates[i] = qexp(runif(1, pexp(0,(rate.val[i])), pexp(1,rate.val[i])), rate.val[i])
  }
  diag(Pol) = selfrates
  
  return(list(Pol=Pol, gen_Load=gen_Load))
}
gen_Success_Matrix = cmpfun(gen_Success_Matrix)