Grid_Set_Match_Runner <- function(numSims, dims,...) {
  
  t1 = proc.time()
  #set up a bunch of empty things for storing data from the simulations.
  field = replicate(numSims, matrix(0, dims, dims), simplify=F)
  gen_Load = replicate(numSims, matrix(0, dims, dims), simplify=F)
  peak_Bloom = replicate(numSims, matrix(0, dims, dims), simplify=F)
  sd_Bloom = replicate(numSims, matrix(0, dims, dims), simplify=F)
  area.occ = replicate(numSims, matrix(0, dims, dims), simplify=F)
  fruit = replicate(numSims, matrix(0, dims, dims), simplify=F)
  abortmat = replicate(numSims, matrix(0, dims, dims), simplify=F)
  probes = replicate(numSims, matrix(0, dims, dims), simplify=F)
  probes.viable = replicate(numSims, matrix(0, dims, dims), simplify=F)
  self.visits = replicate(numSims, matrix(0, dims, dims), simplify=F)
  only.self = replicate(numSims, matrix(0, dims, dims), simplify=F)
  fruitset = replicate(numSims, matrix(0, dims, dims), simplify=F)

  yield = rep(0, numSims)
  globalabort = rep(0, numSims)
  globalfruitset = rep(0, numSims)
  
  for (i in 1:numSims) { #Main loop that calls the simulation and tracks the results from each run.
    tmp = Grid_Set_Match()
    
    field[[i]] = tmp$field
    gen_Load[[i]] = tmp$gen_Load
    peak_Bloom[[i]] = tmp$peak_Bloom
    sd_Bloom[[i]] = tmp$sd_Bloom
    area.occ[[i]] = tmp$area.occ
    fruit[[i]] = tmp$fruit
    abortmat[[i]] = tmp$abortmat
    probes[[i]] = tmp$probes
    probes.viable[[i]] = tmp$probes.viable
    self.visits[[i]] = tmp$self.visits
    only.self[[i]] = tmp$only.self
    fruitset[[i]] = tmp$fruitset
                                                                     
    yield[i] = tmp$yield
    globalabort[i] = tmp$globalabort
    globalfruitset[i] = tmp$globalfruitset
    
    print(sprintf('Simulation %d complete', i))
    
  }
  t2 = proc.time()
  return(list(field=field, gen_Load=gen_Load, peak_Bloom=peak_Bloom, sd_Bloom=sd_Bloom, area.occ=area.occ, fruit=fruit, 
              abortmat = abortmat, probes=probes, self.visits=self.visits, only.self = only.self,
              fruitset = fruitset, yield = yield, globalabort=globalabort, globalfruitset=globalfruitset,
              probes.viable = probes.viable, time=(t2[3]-t1[3])/60))
}
Grid_Set_Match_Runner = cmpfun(Grid_Set_Match_Runner)
Blueberry_Sim_Baseline = Grid_Set_Match_Runner(1000, 15)
save(Blueberry_Sim_Baseline, file="Blueberry_Sim_Baseline.Rdata")