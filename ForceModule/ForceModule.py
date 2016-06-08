import numpy

# Read in the positions of the target atoms, find out total number of atoms
def ExtractPositions(filename,targets):
    
    File = open(filename, "r")
    line_split = File.readlines()
    File.close()
    
    z_list=[]
    for i in targets:
        z = float(line_split[i-1].split()[4])
        z_list.append(z)
        
    return z_list, len(line_split)-3

# Write the forces to file
def WriteForces(filename,forces,targets,numatoms):
    
    targets = [i for i in targets]+[-1] # Need filler at the end for comparison when index matching
    
    File = open(filename, "w")
    TL_index=0; # Target List Index of various lists at which atomid, z coordinate, or force will be found
    for i in numpy.arange(0,numatoms):
        if i+1==targets[TL_index]: # When the index matches an atomid that needs a force
            File.write(`i+1`+' '+'0'+' '+`0.0`+' '+`0.0`+' '+`forces[TL_index]`+'\n') # Write that index and the force
            TL_index+=1
        else:
            File.write(`i+1`+' '+'0'+' '+`0.0`+' '+`0.0`+' '+`0.0`+'\n') # If no force is needed, write the index and zero force.
    File.write('0.0') 
    File.close()
    
# Construct the simulation and experimental densities
def DensCalc(ExpFile,z_list,step=None,zmin=None,zmax=None):
    
    File = open(ExpFile,"r")
    line_split = File.readlines()
    File.close()
    expmin = float(line_split[0].split()[0])
    expmax = float(line_split[0].split()[1])
    expstep = float(line_split[0].split()[2])
    ExpDens = []
    for i in line_split[1:]:
        try:
            ExpDens.append(float(i))
        except:
            pass
    if step == None:
        step = expstep
    if zmin == None:
        zmin = min(expmin,min(z_list))-10*step
    if zmax == None:
        zmax = max(expmax,max(z_list))+10*step

    zrange = numpy.arange(zmin,zmax+step,step)
    
    # Deal with a mismatch of defining grids
    EDens = []
    for z in zrange:
        if z<=expmin or z>=expmax:
            EDens.append(0)
        else:
            indx = int(numpy.floor((z-expmin)/expstep))
            EDens.append(ExpDens[indx]+(ExpDens[indx+1]-ExpDens[indx])*(z-expmin-indx*expstep)/expstep)
    ExpDens=numpy.array(EDens)
    
    SimDens=0*zrange
    for z in z_list:
        SimDens = SimDens+((2*numpy.pi*1.5**2)**-0.5)*numpy.exp(-0.5*(zrange-z)**2./1.5**2.)
                
    return zrange,ExpDens,SimDens
    
# Calculating the Force
def ForceCalc(zrange,ExpDens,SimDens,z_list,func=None):
    if func == None:
        def func(ed,sd):
            return sd-ed
        
    U = func(ExpDens,SimDens)
     
    zstep = zrange[1]-zrange[0]
    forces = []
    for z in z_list:
        indx = int(numpy.floor((z-zrange[0])/zstep))
        forces.append(-1.*(U[indx+1]-U[indx])/zstep)
    
    return forces

if __name__ == "__main__":
    import sys
    
    ARGS=sys.argv[1:]
    targets = eval(ARGS[3])
    z_list,numatoms = ExtractPositions(ARGS[0],targets)
    zrange,ExpDens,SimDens = DensCalc(ARGS[2],z_list)
    forces = ForceCalc(zrange,ExpDens,SimDens,z_list)
    WriteForces(ARGS[1],forces,targets,numatoms)
