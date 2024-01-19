#The gen_Activity_Matrix function creates a matrix that stores the number of turns each pollinator will get during each
#hour of the model, as dictated by that pollinator type's average activity level and the prevaling weather conditions
#during that hour.

gen_Activity_Matrix = function(num_Hours, num_Days, mean_Weather_Qual, sd_Weather_Qual, sd_Weather_Qual_Change,
                               mean_Turns, sd_Turns, num_Agents, PollTypeVec) {
  
  #For each pollinator type, each pollinator gets a random base number of moves per hour that is strongly influenced
  #by parameters specified by the user.
  MovesMatrix = matrix(rep(0, num_Days*num_Agents*num_Hours), num_Hours*num_Days, num_Agents)
  for (bee in 1:num_Agents) {
    mean1 = mean_Turns[PollTypeVec[bee]]
    sd1 = sd_Turns[PollTypeVec[bee]]
    MovesMatrix[,bee] = round(rtruncnorm(num_Hours*num_Days, 0, mean1+(4*sd1), mean1, sd1))
  }
  
  #The function then goes type by type and determines how good or bad for foraging the weather will be on each day of the
  #model run, and then how bad each hour will be within each day. It does this by drawing random quality levels for the first
  #hour of each day and then making the quality of each subsequent hour in each day a modulation of the previous hour's 
  #quality.
  Types = unique(PollTypeVec)
  schedule = matrix(rep(0, num_Agents*num_Hours*num_Days), num_Hours*num_Days, num_Agents)
  for (type in Types) {
    dayqual = rtruncnorm(num_Days, 0, 1, mean_Weather_Qual[type], sd_Weather_Qual[type])
    weathermat = matrix(rep(0, num_Days*num_Hours), num_Hours, num_Days)
    weathermat[1,] = dayqual
    for (day in 1:num_Days) {
      for(hour in 2:num_Hours) {
        last = weathermat[hour-1, day]
        weathermat[hour,day] = rtruncnorm(1, 0, 1, last, sd_Weather_Qual_Change[type])
      }
    }
    weathervec = as.vector(weathermat)
    weathervec = rep(weathervec, each = length(which(PollTypeVec == type)))
    weathermodmat = matrix(weathervec, num_Hours*num_Days, length(which(PollTypeVec == type)), byrow=T)
    col.vals = which(PollTypeVec == type)
    schedule[,col.vals] = round(MovesMatrix[,col.vals] * weathermodmat)
  }
  
  return(list(schedule = schedule, weathermodmat = weathermodmat, MovesMatrix = MovesMatrix))
}
gen_Activity_Matrix = cmpfun(gen_Activity_Matrix)