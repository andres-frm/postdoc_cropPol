#The gen_Bloom_Schedule function creates a unique bloom phenology for every plant in the grid, with the restriction that
#all the flowers the plant will produce will open during the simulation run (i.e. the plant's phenology may be truncated).

gen_Bloom_Schedule = function (num_Days, grid_Dim, mean_Peak_Bloom_Day, sd_Peak_Bloom_Day, mean_Bloom_Range, 
    sd_Bloom_Range, field) {
  
  #The function starts by drawing random values for the peak bloom day for each plant and for the standard deviation around 
  #that peak bloom day, which will affect the "peaked" the bloom phenology ends up being.
  means = matrix(rnorm(grid_Dim*grid_Dim, mean_Peak_Bloom_Day, sd = sd_Peak_Bloom_Day), nrow=grid_Dim)
  sds = matrix(rnorm(grid_Dim*grid_Dim, mean_Bloom_Range, sd = sd_Bloom_Range), nrow=grid_Dim)

  #Next, the function creates several arrays that have a "sheet" for every day being modeled. These arrays spread out the
  #data generated above so that the function can calculate what proportion of a plant's flowers should open on each day. 
  #The last step does this by turning every vector "down the sheets" into a density function with the mean and 
  #standard deviation assigned to that corresponding cell.
  meansArray=array(rep(means,num_Days),c(grid_Dim,grid_Dim,num_Days))
  sdsArray=array(rep(sds,num_Days),c(grid_Dim,grid_Dim,num_Days))
  notallowed = which(sdsArray < 0); sdsArray[notallowed] = 0
  propsArray = array(rep(1:num_Days, each=grid_Dim^2), c(grid_Dim, grid_Dim, num_Days))
  propsArray = array(dnorm(c(propsArray),mean=c(meansArray),sd=c(sdsArray)),c(grid_Dim,grid_Dim,num_Days))
  
  #The density functions above will not sum to one for any vector "down the sheets," so next the function sums each cell's 
  #"down the sheets" vector and divides every proportion therein by the total so that the result sums to one. This is the 
  #step where the truncation mentioned above occurs, as 100% of each plant's flowers are coerced to open during the bloom
  #period being modeled.
  propsArraySumArray=array(rep(apply(propsArray,c(1,2),sum),num_Days),c(grid_Dim,grid_Dim,num_Days)) 
  propsArray = propsArray/propsArraySumArray 
  
  #Lastly, the function takes the proportion of flowers that is supposed to open for every plant on every day and creates a 
  #new array using each plant's flower totals to back-translate the proportions into actual flower totals. It does this by
  #going day by day (i.e. sheet by sheet). Because flower numbers need to be whole numbers, this process can fail to allot, or
  #allot too many, flowers. The function corrects for this by checking the number of flowers alloted to each plant up to this
  #point to the number that "should" have been alloted, and redistributes flowers among the preceding days to correct for
  #any discrepancies.
  cumFlowersArray = array(0, c(grid_Dim,grid_Dim,num_Days)) 
  cumFlowersArray[,,1] = round(propsArray[,,1]*field) 
  flowersArray = array(0, c(grid_Dim,grid_Dim,num_Days)) 
  flowersArray[,,1] = cumFlowersArray[,,1] 
  for (day in 2:num_Days) {  
    cumTotalProp = apply(propsArray[,,1:day],c(1,2),sum) 
    cumFlowersArray[,,day] = round(cumTotalProp*field) 
    flowersArray[,,day] = cumFlowersArray[,,day]-cumFlowersArray[,,day-1] 
  }
  
  return(list(flowersArray=flowersArray, means=means, sds=sds))
}
gen_Bloom_Schedule = cmpfun(gen_Bloom_Schedule)
