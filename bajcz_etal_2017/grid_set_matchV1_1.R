#The simulation model requires a couple of custom packages. These will need to be installed
#on your computer in your working directory and turned on. The code below probably will not work on your computer for that reason.
library("abind", lib.loc="~/R/win-library/3.3")
library("compiler", lib.loc="C:/Program Files/R/R-3.2.4/library")
library("truncnorm", lib.loc="~/R/win-library/3.2")
library("GeneralizedHyperbolic", lib.loc="~/R/win-library/3.3")

#Formally establishes the main simulating function called Grid_Set_Match.
Grid_Set_Match = function(num_Agents = sum(num_Each_Poll_Type),
                          num_Each_Poll_Type = c(36,6,8),
                          mean_Superagents = 402, 
                          sd_Superagents = 0,
                          num_Days=29,
                          num_Hours=12, 
                          grid_Dim = 15, 
                          use_Custom_Grid=FALSE, 
                          custom_Grid = 0, 
                          shape_Flowers_Per_Node = 15.865, 
                          rate_Flowers_Per_Node = 3.19, 
                          shape_Nodes_Per_M2 = 2.81, 
                          rate_Nodes_Per_M2 = 0.00155, 
                          chi_M2_Occupied=9.65, 
                          psi_M2_Occupied=0.0164, 
                          lambda_M2_Occupied=0.033,
                          mean_M2_Occupied = 49.63,
                          flowers_Polls_Ratio = 5012, 
                          stapial_Auto_Prop = 0, 
                          redistribute_Prop = 0,
                          flower_Reduction_Prop=1, 
                          flowers_At_Nest = FALSE,
                          mean_Peak_Bloom_Day = 15, 
                          sd_Peak_Bloom_Day = 4.25, 
                          mean_Bloom_Range = 10.5,
                          sd_Bloom_Range = 2, 
                          num_Attractive_Days = 15, 
                          num_Viable_Days = 7, 
                          finite_Flower_Longevity = TRUE, 
                          mean_Crossing_Prob = 0.418, 
                          sd_Crossing_Prob = 0.154,
                          max_Pollen_Deposited = 25,
                          selfing_Prob_Adjust = 0.106,
                          mean_Turns = c(134, 450, 257), 
                          sd_Turns = c(58, 120, 80), 
                          mean_Weather_Qual = c(0.506, 0.952, 0.732), 
                          sd_Weather_Qual = c(0.1, 0.0428, 0.0886), 
                          sd_Weather_Qual_Change = c(0.130, 0.0317, 0.119),
                          mean_Pollen_Deposited = c(3.5, 26, 19.5),
                          sd_Pollen_Deposited = c(1.5, 5, 4.5), 
                          turns_Per_Bout = c(420, 1944, 240), 
                          nest_Return_Deduction = 0.5, 
                          nest_Locations = list(type1=c(113), 
                                                type2=sample(1:grid_Dim^2, size=num_Each_Poll_Type[2], replace=T), 
                                                type3=sample(1:grid_Dim^2, size=num_Each_Poll_Type[3], replace=T)), 
                          no_Switch_Prob = c(0.833, 0.969, 0.857),  
                          xoffsets1 = c(-1, 0, 1, -1, 0, 1, -1, 0, 1), 
                          yoffsets1 = c(1, 1, 1, 0, 0, 0, -1, -1, -1),
                          xoffsets2 = c(-1, 0, 1, -1, 0, 1, -1, 0, 1), 
                          yoffsets2 = c(1, 1, 1, 0, 0, 0, -1, -1, -1),
                          xoffsets3 = c(-1, 0, 1, -1, 0, 1, -1, 0, 1), 
                          yoffsets3 = c(1, 1, 1, 0, 0, 0, -1, -1, -1), 
                          max_Carryover = 25, 
                          intrahive_Transfer = TRUE, 
                          max_Pollen_Transfer = 5,
                          plotting_Freq=0) { #Plotting is currently disabled
  
  #Starts the timer so that the duration of each simulation can be determined.
  t1 = proc.time()
  
  #Unpacks type information about each pollinator being modeled so the function can handle it more easily.
  num_Type1 = num_Each_Poll_Type[1]
  num_Type2 = num_Each_Poll_Type[2]
  num_Type3 = num_Each_Poll_Type[3]
  
  #Immediately checks for several common errors that the simulation function is not equipped to handle.
  if ((num_Type1 + num_Type2 + num_Type3) != num_Agents) { 
    stop('Illegal error #3! Please make sure poll numbers sum to num_Agents!')
  }
  if (num_Viable_Days > num_Attractive_Days && finite_Flower_Longevity == TRUE) {
    stop('Illegal error #4! Viability cannot exceed num_Attractive_Days!')
  }
  if (use_Custom_Grid == TRUE && custom_Grid == 0) {
    stop('Illegal error #5! Chosen field must be specified if use_Custom_Grid = TRUE!')
  }
  if (use_Custom_Grid == TRUE && length(custom_Grid[,1]) != length(custom_Grid[1,])) {
    stop('Illegal error #6! Chosen field must be a square matrix!')
  }
  
  #If flowers will not lose viability, then max viability is set equal to the number of days being modeled. 
  if (finite_Flower_Longevity == FALSE) {
    num_Attractive_Days = num_Days
    num_Viable_Days = num_Days
  }
  
  #Builds a wrapper of sample that will behave itself better during internal simulation operations when length of the vector
  #being sampled from may occasionally be 1. 
  sample.vec <- function(x, ...) {
    x[sample(length(x), ...)]
  }
  sample.vec = cmpfun(sample.vec)
  
  #If a custom grid is not being used, the simulation calls the gen_Grid function to make a grid for this run.
  defaultflowernum = ((mean_Superagents * num_Agents * flowers_Polls_Ratio)/grid_Dim^2)
  if (use_Custom_Grid == FALSE) {
    setupfield = gen_Grid(defaultflowernum, grid_Dim, stapial_Auto_Prop, shape_Flowers_Per_Node, rate_Flowers_Per_Node, 
        lambda_M2_Occupied, shape_Nodes_Per_M2, rate_Nodes_Per_M2, mean_M2_Occupied, chi_M2_Occupied, psi_M2_Occupied)
    field = setupfield$fieldflr
    area.occ = setupfield$area.occ #Area occupied by each clone for reporting.
    
    #If a custom grid is being used, though, the simulation needs to know a few metrics about that new field.
  } else  {
    field = custom_Grid
    grid_Dim = length(field[,1])
    area.occ = NA #We wouldn't know what this was.
  }
  
  #Error check to make sure that nest locations were not specified illogically.
  if (any(unlist(nest_Locations)>length(field)) == TRUE) {
    stop('Illegal error #7! You cannot specify nest locations that are not in the field!')
  }
  
  #If flowers can't occupy the same squares as the nests of type1 pollinators, as might be the case with honeybees, the
  #model must erase any flowers it put in those locations.
  if (flowers_At_Nest == FALSE) {
    field[nest_Locations[[1]]] = 0
  }
  
  #The simulation calls the gen_Polls_DF function to create a data frame that will store important info about each pollinator.
  setupbees1 = gen_Polls_DF(field, num_Agents, num_Type1, num_Type2, no_Switch_Prob, mean_Pollen_Deposited, 
                            sd_Pollen_Deposited, turns_Per_Bout, mean_Superagents, sd_Superagents, nest_Locations)
  polls = setupbees1$polls
  
  #The simulation calls the gen_Activity_Matrix funtion to create a matrix of turns for each pollinator for every hour
  #of the simulation run.
  weather = gen_Activity_Matrix(num_Hours, num_Days, mean_Weather_Qual, sd_Weather_Qual, sd_Weather_Qual_Change, mean_Turns, 
                                sd_Turns, num_Agents, PollTypeVec = setupbees1$PollType)
  schedule = weather$schedule; prior.schedule = schedule
  
  #Next, we must check to make sure that users aren't trying to further manipulate custom fields in unintended ways.
  if (redistribute_Prop > 0 & use_Custom_Grid == TRUE) { 
    stop("Do not use redistribute_Prop if using a custom field; it will add flowers to spots that should be empty.")
  }
  
  #If any form of homogenization in flowers numbers is being modeled, the next few lines carry these processes out.
  field2 = round(field * redistribute_Prop); allflowers = sum(field2); flowersubsidy = round(allflowers/(grid_Dim^2))
  field = field + flowersubsidy; field = field - field2
  field = round(field * flower_Reduction_Prop) 
  
  #Next, the model creates tons of empty vectors/matrices to store data in as it is being gathered.
  fruit = matrix(0, grid_Dim, grid_Dim) 
  abortmat = matrix(0, grid_Dim, grid_Dim) 
  fruitset = matrix(0, grid_Dim, grid_Dim) 
  flowers.Move = matrix(0, grid_Dim, grid_Dim) 
  flowers.Pol = matrix(0, grid_Dim, grid_Dim) 
  PollMoveCounter = rep(0, num_Agents) 
  probes = matrix(0, grid_Dim, grid_Dim) 
  probes.viable = matrix(0, grid_Dim, grid_Dim) 
  self.visits = matrix(0, grid_Dim, grid_Dim) 
  only.self = matrix(0, grid_Dim, grid_Dim) 
  
  #Similarly, it now sets up baseline values for various counters that will track the progression of certain processes like the
  #passage of time.
  #displayCounter = 0 #Deprecated feature
  HourCounter = 0 
  DayCounter = 1 
  
  #Every plant needs its own unique ID number, which is given next.
  plantID = matrix(seq(1:grid_Dim^2), grid_Dim, grid_Dim) 
  
  #Each pollinator's pollen history is tracked in a matrix called pollen, which is created next.
  if (intrahive_Transfer == TRUE & max_Pollen_Transfer > max_Carryover) {
    pollen.length = max_Pollen_Transfer
  } else { pollen.length = max_Carryover}
  pollen = matrix(rep(0, num_Agents * pollen.length), num_Agents, pollen.length) 
  
  #Each pollinator type has its own "movement neighborhood" that determines which grid cells can be reached during movement.
  xlist = list(xoffsets1, xoffsets2, xoffsets3)
  ylist = list(yoffsets1, yoffsets2, yoffsets3)
  
  #The simulation next calls the function gen_Bloom_Schedule to create a unique phenology of bloom for each plant in the grid.
  flowersArray1 = gen_Bloom_Schedule(num_Days, grid_Dim, mean_Peak_Bloom_Day, sd_Peak_Bloom_Day, mean_Bloom_Range, sd_Bloom_Range, field)
  flowersArray = flowersArray1$flowersArray
  peak_Bloom = flowersArray1$means; sd_Bloom = flowersArray1$sds #Interesting data for reporting.
  
  #Pollinators will only be able to interact with flowers that are open at that point, so now we "open" the first day's flowers.
  currentflowersArray = array(flowersArray[,,1], dim = c(grid_Dim,grid_Dim,DayCounter))
  
  #The simulations now calls the function gen_Success_Matrix to create a probability of pollination success for every possible
  #pollen donor/receiver pairing. 
  successes = gen_Success_Matrix(mean_Crossing_Prob, sd_Crossing_Prob, grid_Dim, selfing_Prob_Adjust)
  Pol = successes$Pol
  gen_Load = successes$gen_Load
  
  #The simulation now begins in earnest. Most externally, the model goes hour by hour to simulate events. 
  for (hour in 1:length(schedule[,1])) { 
    HourCounter = HourCounter + 1 
    
    #However, days are also relevant time steps for certain processes, so the simulation tracks those as well.
    if (DayCounter == 1 && hour == 1) {
      print(sprintf('Starting day %g', DayCounter))
    }
    
    #If we've entered a new day, the simulation performs several updating tasks specific to the change of days.
    if (HourCounter == (num_Hours + 1)) { 
      DayCounter = DayCounter + 1 
      polls = setupbees1$polls 
      HourCounter = 1 
      PollMoveCounter = PollMoveCounter * 0 
      pollen = pollen * 0 #reset the pollen matrix
      print(sprintf('Starting day %g', DayCounter))
      
      #The most important of these is updating the matrix of flowers that pollinators can interact with.
      if (DayCounter > num_Attractive_Days && num_Attractive_Days > 1) {
        currentflowersArray = abind(currentflowersArray[,,2:num_Attractive_Days], flowersArray[,,DayCounter], along=3)
      } else {
        if (num_Attractive_Days !=1) {
          currentflowersArray = abind(currentflowersArray, flowersArray[,,DayCounter], along=3)
        } else {
          currentflowersArray = array(flowersArray[,,DayCounter], dim=c(grid_Dim, grid_Dim))
        }
      }
    } 
    
    #Because flowers can be open but non-viable, the simulation uses two seperate matrices to track flowers of both
    #types so that movement and pollination can use one or the other.
    if (is.na(dim(currentflowersArray)[3]) == TRUE) {
      flowers.Move = currentflowersArray
      flowers.Pol = currentflowersArray
    } else {
      flowers.Move = apply(currentflowersArray, c(1,2), sum, na.rm=T)
      flowers.Pol = apply(currentflowersArray[,,seq(to = dim(currentflowersArray)[3], 
          length = min(dim(currentflowersArray)[3], num_Viable_Days))], c(1,2), sum, na.rm=T)
    }
    
    #The simulation also tracks how many total turns and total viable flowers are still present in the field because if either
    #runs out, the model should automatically advance to the next time step to increase model runtime efficiency.
    totalnumFlowers = sum(flowers.Pol)
    totalTurns = sum(schedule[hour,]) 
    
    #Within each hour, pollinators are allowed to act on the grid until either their turns or the flowers are depleted.
    while (totalnumFlowers > 0 && totalTurns > 0) { 
      
      #Pollinators are chosen to act in a random order each time to as closely simulate "simultaneous" action as possible.
      for (poll in sample.vec(which(schedule[hour,]>0))) { 
        if(totalnumFlowers > 0) {
          
          #Once a pollinator's bout is over, the simulation returns it to its nest and performs other update actions.
          #These may include returning the pollinator to its home, reseting its bout counter, subtracting some future turns 
          #from its activity schedule, reseting its pollen history, and adding new pollen to its history.
          if (PollMoveCounter[poll] >= polls$VisitsCounter[poll]) {
            findhome = arrayInd(polls$homes[poll], dim(field))
            polls$x[poll] = findhome[1]; polls$y[poll] = findhome[2]
            PollMoveCounter[poll] = 0
            moves.to.deduct = round(mean_Turns[polls$type[poll]]*nest_Return_Deduction)
            amt.to.deduct = min(schedule[hour,poll], moves.to.deduct)
            schedule[hour,poll] = schedule[hour,poll] - amt.to.deduct
            moves.to.deduct = moves.to.deduct - amt.to.deduct
            timer = 1
            while (moves.to.deduct > 0) {
              if ((hour + timer) > num_Hours) { break }
              amt.to.deduct = min(schedule[hour+timer,poll], moves.to.deduct)
              schedule[hour+timer,poll] = schedule[hour+timer,poll] - amt.to.deduct
              moves.to.deduct = moves.to.deduct - amt.to.deduct
              timer = timer + 1
            }
            if (intrahive_Transfer == TRUE) { 
              which.columns = max(1, (max_Carryover-max_Pollen_Transfer+1))
              this.type = polls$type[poll] 
              all.this.type = which(polls$type == this.type)
              possible.clones = pollen[which(pollen[all.this.type,] > 0)]
              pollen[poll,which.columns:length(pollen[1,])] = sample.vec(possible.clones, max_Pollen_Transfer, replace = T) 
              if (max_Carryover > max_Pollen_Transfer) {
                pollen[poll, 1:max(which.columns-1, 1)] = pollen[poll, 1:max(which.columns-1, 1)] * 0 
              }
            } else {
              pollen[poll,] = pollen[poll,] * 0
            } 
          } 
          
          #Because taking away turns from pollinators affects the duration of the while loop, the simulation checks to make
          #sure the model advances to the next time step appropriately.
          totalTurns = sum(schedule[hour,]) 
          if (schedule[hour,poll] == 0) { next }
          
          #The simulation now enters the movement phase of the current turn. Pollinators are either moved to a new grid cell
          #or remain in the current grid cell. The next block of code accomplishes this. First, a draw is performed to 
          #determine if the pollinator switches grid cells. If it does, the model figures out what grid cells it can reach,
          #and it determines which of these have flowers and how many flowers each has. It uses that information to move
          #the pollinator semi-randomly to a new grid cell it can reach, with cells with more flowers being more likely to
          #be chosen. A few error checks that should NEVER trigger are also in here. This part concludes by identifying
          #which cell the pollinator is now in.
          if ((flowers.Move[polls$x[poll],polls$y[poll]] == 0) || (runif(1) >= polls$stayprob[poll])) {
            newxs = polls$x[poll] + xlist[[polls$type[poll]]] 
            newys = polls$y[poll] + ylist[[polls$type[poll]]]
            coords = cbind(x=newxs, y=newys) 
            keep = which(coords[,1]>=1 & coords[,1]<=grid_Dim & coords[,2]>=1 & coords[,2]<=grid_Dim)
            flrs = flowers.Move[coords[keep,]]
            totalflr = sum(flrs) 
            if (totalflr > 0) {
              probflr = flrs/totalflr
              if (any(flrs < 0) == TRUE) {
                stop('Illegal error #8! Somehow, one or more plants has a negative number of flowers!')
              }
              movedir = sample(length(flrs), 1, prob=probflr) 
              movedir = keep[movedir]
              polls$x[poll] = polls$x[poll] + xlist[[polls$type[poll]]][movedir]
              polls$y[poll] = polls$y[poll] + ylist[[polls$type[poll]]][movedir]
            }  else { 
              tmp1 = which(flowers.Move > 0) 
              if (length(tmp1) > 0) { 
                tmp2 = sample.vec(tmp1, 1)  
                Next = arrayInd(tmp2, dim(flowers.Move)) 
                polls$x[poll] = Next[1] 
                polls$y[poll] = Next[2] 
              } else { 
                stop('Illegal error #1!  Tried to move with no flowers in entire field!')
              } 
            } 
          } 
          newID = plantID[polls$x[poll],polls$y[poll]]
          
          #Now, the simulation enters the pollination phase. First, it checks to see what pollen the pollinator is carrying.
          #Provided that it has some, it determines how many pollen grains it will deposit on the current flower, and it
          #identifies the donors of those pollen grains.
          subset = which(pollen[poll,] > 0) 
          options = length(subset) 
          if (options > 0) { 
            PollenGrains = round(rtruncnorm(1, 0, polls$meanpollen[poll]*4, polls$meanpollen[poll], polls$sdpollen[poll]))
            if (PollenGrains > 0)  {
              probs = (seq(1:options))/sum(seq(1:options))
              pick = sample(1:options, PollenGrains, prob=probs, replace = TRUE)
              pick = subset[pick]
              oldIDs = pollen[poll, pick] 
              
              #The simulation tracks how many probes each plant receives as well as how many probes that result in 
              #self-pollination have occured. 
              probes[newID] = probes[newID] + (1*polls$numIndiv[poll])
              self.pollen = which(oldIDs == newID)
              if (length(self.pollen) > 0) { 
                self.visits[newID] = self.visits[newID] + (1*polls$numIndiv[poll])
              }
              if (length(self.pollen) == length(oldIDs)) { 
                only.self[newID] = only.self[newID] + (1*polls$numIndiv[poll])
              }
              if (flowers.Pol[polls$x[poll], polls$y[poll]] > 0) {
                probes.viable[newID] = probes.viable[newID] + (polls$numIndiv[poll])
                
                #Next, the simulation calculates how likely this particular pollination event is to be successful.
                SuccessProb = mean(Pol[oldIDs, newID]) * (PollenGrains/max_Pollen_Deposited)
                
                #A draw occurs to see if the event was successful. If it was, flowers need to be deducted from the current
                #plant's pool. This is a complicated process. The model first figures out how many flowers are currently
                #viable, then determines how many individuals the current pollinator represents, and it deducts up to that
                #many flowers from the current plant. 
                if (runif(1) < SuccessProb) {
                  if (is.na(dim(currentflowersArray)[3]) == FALSE) { 
                    howmany2 = which(currentflowersArray[polls$x[poll], polls$y[poll], seq(to = dim(currentflowersArray)[3], 
                        from = max(1, dim(currentflowersArray)[3]-num_Viable_Days+1), by=1)] > 0)
                    these.sheets = howmany2 + max(0, dim(currentflowersArray)[3]-num_Viable_Days)
                    currentflowers = currentflowersArray[polls$x[poll], polls$y[poll], these.sheets]
                    sumcurrentflowers = sum(currentflowers)
                    flowers.todo = min(sumcurrentflowers, polls$numIndiv[poll])
                    this.many.flowers = flowers.todo
                    while (flowers.todo > 0) { 
                      howmany2 = which(currentflowersArray[polls$x[poll], polls$y[poll], seq(to = dim(currentflowersArray)[3], 
                          from = max(1, dim(currentflowersArray)[3]-num_Viable_Days+1), by=1)] > 0)
                      these.sheets = howmany2 + max(0, dim(currentflowersArray)[3]-num_Viable_Days) 
                      currentflowers = currentflowersArray[polls$x[poll], polls$y[poll], these.sheets]
                      splits1 = rmultinom(n=1, flowers.todo, prob=currentflowers/sum(currentflowers)) 
                      toomany = which(splits1 > currentflowers)
                      if (length(toomany) >= 1) { 
                        splits1[toomany] = currentflowers[toomany]
                      } 
                      currentflowersArray[polls$x[poll], polls$y[poll],these.sheets] =
                        currentflowersArray[polls$x[poll], polls$y[poll],these.sheets] - splits1
                      flowers.todo = flowers.todo - sum(splits1)
                    }
                  } else {
                    currentflowers = currentflowersArray[polls$x[poll], polls$y[poll]]
                    currentflowersArray[polls$x[poll], polls$y[poll]] = currentflowersArray[polls$x[poll], 
                        polls$y[poll]] - min(polls$numIndiv[poll], currentflowers)
                  }
                  
                  #Once the flowers have been removed from the current plant's pool, the various other matrices that
                  #are tracking the current pool of flowers with which pollinators interact need to be updated. In addition,
                  #fruits need to be added to the current plant's pool.
                  flowers.Move[polls$x[poll],polls$y[poll]] = flowers.Move[polls$x[poll], polls$y[poll]] - this.many.flowers
                  flowers.Pol[polls$x[poll],polls$y[poll]] = flowers.Pol[polls$x[poll], polls$y[poll]] - this.many.flowers
                  totalnumFlowers = totalnumFlowers - this.many.flowers
                  fruit[polls$x[poll],polls$y[poll]] = fruit[polls$x[poll],polls$y[poll]] + this.many.flowers  
                  
                #However, if the pollination event was not successful, the simulation draws to see if the current flower
                #will have been permanently ruined by that incompatible pollen. If it is ruined, the process proceeds 
                #identically to above, except that no fruits are produced at the end. 
                } else {
                  if (PollenGrains > (max_Pollen_Deposited - 1)) { 
                    abortprob = 1
                  } else {
                    abortprob = (PollenGrains/max_Pollen_Deposited)
                  }
                  if (runif(1) < abortprob) {
                    if (is.na(dim(currentflowersArray)[3]) == FALSE) {
                      howmany2 = which(currentflowersArray[polls$x[poll], polls$y[poll], seq(to = dim(currentflowersArray)[3], 
                          from = max(1, dim(currentflowersArray)[3]-num_Viable_Days+1), by=1)] > 0)
                      these.sheets = howmany2 + max(0, dim(currentflowersArray)[3]-num_Viable_Days)
                      currentflowers = currentflowersArray[polls$x[poll], polls$y[poll], these.sheets]
                      sumcurrentflowers = sum(currentflowers)
                      flowers.todo = min(sumcurrentflowers, polls$numIndiv[poll])
                      this.many.flowers = flowers.todo 
                      while (flowers.todo > 0) {
                        howmany2 = which(currentflowersArray[polls$x[poll], polls$y[poll], seq(to = dim(currentflowersArray)[3], 
                            from = max(1, dim(currentflowersArray)[3]-num_Viable_Days+1), by=1)] > 0)
                        these.sheets = howmany2 + max(0, dim(currentflowersArray)[3]-num_Viable_Days) 
                        currentflowers = currentflowersArray[polls$x[poll], polls$y[poll], these.sheets]
                        splits1 = rmultinom(n=1, flowers.todo, prob=currentflowers/sum(currentflowers)) 
                        toomany = which(splits1 > currentflowers)
                        if (length(toomany) >= 1) {
                          splits1[toomany] = currentflowers[toomany]
                        }
                        currentflowersArray[polls$x[poll], polls$y[poll],these.sheets] =
                          currentflowersArray[polls$x[poll], polls$y[poll],these.sheets] - splits1
                        flowers.todo = flowers.todo - sum(splits1)
                      }
                      
                    } else {
                      currentflowers = currentflowersArray[polls$x[poll], polls$y[poll]]
                      currentflowersArray[polls$x[poll], polls$y[poll]] = currentflowersArray[polls$x[poll], 
                        polls$y[poll]] - min(1*polls$numIndiv[poll], currentflowers)
                    }
                    flowers.Move[polls$x[poll],polls$y[poll]] = flowers.Move[polls$x[poll], polls$y[poll]] - this.many.flowers
                    flowers.Pol[polls$x[poll],polls$y[poll]] = flowers.Pol[polls$x[poll], polls$y[poll]] - this.many.flowers
                    totalnumFlowers = totalnumFlowers - this.many.flowers
                    abortmat[polls$x[poll], polls$y[poll]] = abortmat[polls$x[poll], polls$y[poll]] + this.many.flowers
                  }
                } 
              }
            }
          }
          
          #Now, the simulation enters a brief updating phase. Here, the pollinator's bout counter is updated, and its pollen
          #history is updated as well.
          PollMoveCounter[poll] = PollMoveCounter[poll] + 1
          if (max_Carryover >= max_Pollen_Transfer) {
            if (max_Carryover == 1) {
              pollen[poll,1] = newID 
            } else {
              pollen[poll,] = c(pollen[poll,2:length(pollen[1,])], newID)
            }
            } else {
              if (max_Carryover == 1) {
                pollen[poll,] = c(pollen[poll, 1:(length(pollen[1])-1)], newID)
              } else {
              tmp5 = max_Pollen_Transfer - max_Carryover
              first.chunk = pollen[poll,1:tmp5]
              tmp6 = tmp5 + 2
              second.chunk = pollen[poll,tmp6:max_Pollen_Transfer]
              pollen[poll,] = c(first.chunk, second.chunk, newID)
              }
            }
          
          #A levelplot of the current fruit set matrix is displayed after a number of turns specified by the user.
          #This feature is deprecated.
          #displayCounter = displayCounter + 1
          #fruitset = fruit/field
          #blah = which(fruitset=='NaN')
          #fruitset[blah] = NA
          #if(displayCounter==plotting_Freq) {
          #print(levelplot(fruitset))
          #displayCounter=0
          #}
          
        }
      }
      
      #Lastly, the number of turns each pollinator has left, as well as the total number of turns left this hour for all
      #pollinators, are updated to keep the while loop operating.
      totalTurns = totalTurns - length(which(schedule[hour,]>0))
      schedule[hour,which(schedule[hour,]>0)] = schedule[hour,which(schedule[hour,]>0)] - 1 
      
    }
  }
  
  #Before the simulation model ends, it calculates a couple of global data metrics based on the data collected during the run.
  fruitset = fruit/field
  meanfruitset = mean(fruitset, na.rm=TRUE)
  globalfruitset = sum(fruit)/sum(field)
  abortprop = abortmat/field
  globalabort = sum(abortmat)/sum(field)
  yield = sum(fruit)
  
  #It closes by recording the end time so that total simulation run time can be calculated.
  t2 = proc.time()
  return(list(field=field, meanfruitset = meanfruitset, yield=yield, fruit=fruit, fruitset=fruitset, 
              time=((t2[3]-t1[3])/60), abortmat=abortmat, abortprop=abortprop, Pol = Pol, self.visits = self.visits,
              pollen=pollen, schedule=prior.schedule, globalfruitset = globalfruitset, probes = probes,
              globalabort=globalabort, only.self = only.self, probes.viable = probes.viable, area.occ=area.occ, gen_Load=gen_Load,
              peak_Bloom=peak_Bloom, sd_Bloom = sd_Bloom))
}
Grid_Set_Match=cmpfun(Grid_Set_Match)