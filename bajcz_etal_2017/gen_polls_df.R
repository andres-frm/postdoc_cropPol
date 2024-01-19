#The gen_Polls_DF function creates a data frame that will store important identifying and type information about each 
#pollinator being modeled. The function exists moreso to consolate disparate information and has relatively few calculations
#to perform.

gen_Polls_DF = function (field, num_Agents, num_Type1, num_Type2, no_Switch_Prob, mean_Pollen_Deposited, sd_Pollen_Deposited, 
                         turns_Per_Bout, mean_Superagents, sd_Superagents, nest_Locations) {
  
  #First, the function determines how many "real" individuals each pollinator agent will represent.
  numIndiv = round(rtruncnorm(num_Agents, 0, Inf, mean_Superagents, sd_Superagents))
  
  #Then, several empty vectors are created to store functional information needed to properly align information by type.
  stayProb = rep(0, num_Agents)
  meanPollen = rep(0, num_Agents)
  sdPollen = rep(0, num_Agents)
  VisitsCounter = rep(0, num_Agents)
  HiveHome = rep(0, num_Agents)
  
  #Pollinators are then assigned to their respective types so that agents of different types can be modeled differently.
  PollType = seq(1:num_Agents)
  type1 = which(PollType > 0 & PollType < (num_Type1 + 1))
  type2 = which(PollType > (num_Type1) & PollType < (num_Type1 + num_Type2 + 1))
  type3 = which(PollType > (num_Type1 + num_Type2) & PollType < (num_Agents + 1))
  PollType[type1] = 1; PollType[type2] = 2; PollType[type3] = 3
  
  #Characteristics that are pollinator-type-specific are then assigned to each pollinator.
  stayProb[type1] = no_Switch_Prob[1]; no_Switch_Prob[type2] = no_Switch_Prob[2]; stayProb[type3] = no_Switch_Prob[3]
  meanPollen[type1] = mean_Pollen_Deposited[1]; meanPollen[type2] = mean_Pollen_Deposited[2]; meanPollen[type3] = mean_Pollen_Deposited[3]
  sdPollen[type1] = sd_Pollen_Deposited[1]; sdPollen[type2] = sd_Pollen_Deposited[2]; sdPollen[type3] = sd_Pollen_Deposited[3]
  VisitsCounter[type1] = turns_Per_Bout[1]; VisitsCounter[type2] = turns_Per_Bout[2]; VisitsCounter[type3] = turns_Per_Bout[3]
  
  #This includes an assignment to a specific nest location. Each pollinator type is allowed to have any number of seperate
  #nest locations within the grid, and the function next assigns each pollinator of each type to one of the nest locations
  #specified for that type.
  for (polltype in unique(PollType)) { 
    thesepolls = which(PollType == polltype) 
    numHives = length(nest_Locations[[polltype]])
    hiveplacement = rep_len(seq(from=1, to=numHives, by=1), length=length(thesepolls))
    Hivespotvec = nest_Locations[[polltype]][hiveplacement]
    HiveHome[thesepolls] = Hivespotvec
  }
  xv = arrayInd(HiveHome, dim(field))[,1]
  yv = arrayInd(HiveHome, dim(field))[,2]
    
  #Lastly, the polls data frame is formally assembled to be passed along to the main simulation model.
  polls = data.frame(x = xv, y = yv, type = PollType, stayprob = stayProb, meanpollen = meanPollen, sdpollen = sdPollen,
       VisitsCounter = VisitsCounter, numIndiv = numIndiv, homes=HiveHome)
  return(list(polls = polls, PollType=PollType))
  
}
gen_Polls_DF = cmpfun(gen_Polls_DF)