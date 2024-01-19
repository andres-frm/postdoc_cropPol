#The gen_Grid function creates a matrix of cells, each associated with a number that represents the total number of flowers
#the plant in that cell will produce over the duration of the simulation. 

gen_Grid = function(defaultflowernum, grid_Dim, stapial_Auto_Prop, shape_Flowers_Per_Node, rate_Flowers_Per_Node, 
    lambda_M2_Occupied, shape_Nodes_Per_M2, rate_Nodes_Per_M2, mean_M2_Occupied, chi_M2_Occupied, psi_M2_Occupied, ...) {
    
    #The function starts by creating a fully homogeneous field, and an empty matrix of the same size to store modulating
    #factors that will eventually be used to adjust the flower totals in the first matrix.
    fieldflr = matrix(defaultflowernum, grid_Dim, grid_Dim)
    fieldgenenv = matrix(0, grid_Dim, grid_Dim)
    
    #The number of flowers each plant is expected to produce is a function of three factors: flowers per node, nodes per area
    #and area occupied, which each follow a specific distribution. Values are drawn for each factor for each plant, and then
    #these values are mean-adjusted to turn them into factors relating how fecund each plant is relative to the average plant.
    #These three factors are multiplied together to get a single factor that relates the fecundity of a clone to average
    #fecundity, and it is this factor that is eventually multiplied by the average flower number. For example, if a plant 
    #is only half as fecund as the average plant for all three factors, it's final factor will be 0.5*0.5*0.5=0.125, which
    #means this plant will produce only 1/8 as many flowers as the average plant. 
    area.occ = matrix(rep(0, grid_Dim^2), grid_Dim, grid_Dim) #Empty matrix to store area data in.
    for (row in 1:grid_Dim) {
      for (column in 1:grid_Dim) {
      flr = 0; node = 0; size = 0
      flr = rgamma(1, shape=shape_Flowers_Per_Node, rate=rate_Flowers_Per_Node)
      flrfac = flr/(shape_Flowers_Per_Node/rate_Flowers_Per_Node)
      node = rgamma(1, shape=shape_Nodes_Per_M2, rate=rate_Nodes_Per_M2)
      nodefac = node/(shape_Nodes_Per_M2/rate_Nodes_Per_M2)
      size = rgig(1, chi=chi_M2_Occupied, psi=psi_M2_Occupied, lambda=lambda_M2_Occupied)
      area.occ[row,column] = size
      sizefac = size/mean_M2_Occupied
      fieldgenenv[row,column] = flrfac * nodefac * sizefac 
      }
    }

#The spatdep function is used to put some spatial autocorrelation into these final factors so that plants closer to each
#other can be modeled as more similar than plants further apart. It does this by spatially averaging the factors to some 
#degree specified by the user.
spatdep = function(A, stapial_Auto_Prop = 0) {
  Lx = dim(A)[1]
  Ly = dim(A)[2]
  north = A[,((1:Ly)+1-1)%%Ly + 1]
  northwest = A[((1:Lx)-1-1)%%Lx + 1,((1:Ly)+1-1)%%Ly + 1]
  south = A[,((1:Ly)-1-1)%%Ly + 1]
  southwest = A[((1:Lx)-1-1)%%Lx + 1,((1:Ly)-1-1)%%Ly + 1]
  east = A[((1:Lx)+1-1)%%Lx + 1,]
  northeast = A[((1:Lx)+1-1)%%Lx + 1,((1:Ly)+1-1)%%Ly + 1]
  west = A[((1:Lx)-1-1)%%Lx + 1,]
  southeast = A[((1:Lx)+1-1)%%Lx + 1,((1:Ly)-1-1)%%Ly + 1]
  return((1-stapial_Auto_Prop)*A + stapial_Auto_Prop/(8*(north+northwest+south+southwest
                                   +east+northeast+west+southeast)))
}
fastspatial = cmpfun(spatdep)

#Any spatial autocorrelation desired is applied, and then the default flower numbers are modulated by the final factors
#calculated above to create a single matrix of flower totals.
    fieldgenenv = fastspatial(fieldgenenv, stapial_Auto_Prop)
    fieldflr = round(fieldflr * fieldgenenv)
  return(list(fieldflr=fieldflr, fieldgenenv=fieldgenenv, area.occ=area.occ))
}
gen_Grid = cmpfun(gen_Grid)
