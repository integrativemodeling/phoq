import dl
import sys
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | dl.RTLD_GLOBAL)

import IMP
import IMP.algebra
import IMP.atom
import IMP.membrane
import IMP.container
import IMP.core
import IMP.rmf
import RMF
import math
import IMP.isd2

# PARAMETERS
W_BIAS_STRIDE=10000
# max step for Monte Carlo moves
MAXOMEGA_=0.2
MAXTRANS_=0.15
MAXROT_=0.03
MAXEPSILON_=0.01
MAXWEIGHT_=0.1

NITER_=700000
NOPT_=500
HOMOLOGY_CUTOFF_=8.0
ELASTIC_CUTOFF_=6.0
W_STRIDE=50
KAPPA_=100.0
NGRID_=1000
# CROSSLINK PARA
# typical beta derived from average standard deviation
# of experimental error
BETA_=0.025
SIGMAG_=0.02
# REM PARAMETERS
TEMPMIN_=0.8
TEMPMAX_=1.2
# WTE PARAMETERS
EMIN_=-4000.0
EMAX_=0.0
SIGMA_=10.0
GAMMA_=16.0
W0_=0.01
# NUMBER OF COPIES
NCOPIES_=2
# MIXING PARAMETER FOR CA FF
NATIVENESS_=0.3
KT_CAFF_=1.0

# MODEL
m = IMP.Model()

# set up replica exchange
rem=IMP.membrane.ReplicaExchange()
# get number of replicas
nproc=rem.get_number_of_replicas()
# create array of temperature
temp_array=rem.create_temperatures(TEMPMIN_,TEMPMAX_,nproc)
# get my replica index      
myindex=rem.get_my_index()

# set initial value for temperature 
tempval=temp_array[myindex]

# set parameters
rem.set_my_parameter("temp",[tempval])

def read_alignments(filename,model,temp):
    f=open(filename,"r")
    align={}
    for line in f.readlines():
        riga=(line.strip()).split()

        try:
         s1=IMP.atom.Selection(model, chains=riga[1], residue_index=int(riga[2]), residue_type=IMP.atom.get_residue_type(riga[0]), atom_type=IMP.atom.AT_CA)
         p1=s1.get_selected_particles()[0]
        except:
         print "residue %d of chain %s not found in model" % (int(riga[2]),riga[1])
         continue
        try:
         s2=IMP.atom.Selection(temp, chains=riga[4], residue_index=int(riga[5]), residue_type=IMP.atom.get_residue_type(riga[3]), atom_type=IMP.atom.AT_CA) 
         p2=s2.get_selected_particles()[0]
        except:
         print "residue %d of chain %s not found in template" % (int(riga[5]), riga[4])
         continue

        align[p1]=p2

    f.close()
    return align

def load_pdb(name,filename):
    prot=IMP.atom.read_pdb(filename,m,IMP.atom.CAlphaPDBSelector()) 
    prot.set_name(name)
    return prot

def create_rigid_body(chain,res0,res1):
    s=IMP.atom.Selection(chain, residue_indexes=[(res0,res1+1)])
    atoms=[]
    for p in s.get_selected_particles():
        atoms.append(IMP.core.XYZR(p))
    prb=IMP.Particle(m)
    rb=IMP.core.RigidBody.setup_particle(prb,atoms) 
    return rb

def get_excluded_volume(prot,kappa,excluded_residues=[]):
    rset=IMP.RestraintSet('Excluded_Volume')
    atoms_to_use=[]
    for atom in IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE):
        residue=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent())
        resindex=residue.get_index()
        if(resindex in excluded_residues): continue
        restype=residue.get_residue_type()
        vol=IMP.atom.get_volume_from_residue_type(restype)
        radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
        IMP.core.XYZR(atom).set_radius(radius)
        atoms_to_use.append(atom)
    lsc=IMP.container.ListSingletonContainer(m)
    lsc.add_particles(atoms_to_use)
    evr=IMP.core.ExcludedVolumeRestraint(lsc,kappa)
    evr.set_name('Excluded_Volume')
    rset.add_restraint(evr)
    return rset

def get_distance_restraint(p0, p1, d0, kappa):
    h=IMP.core.Harmonic(d0,kappa)
    dps=IMP.core.DistancePairScore(h)
    pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p0,p1))
    return pr

def read_potential_dihedral(filename,string,mix=False):
# read potentials for dihedral
   score_dih=[]
   phi0=[]; phi1=[]
   for i in range(0,36):
       phi0.append(i*10.0/180.0*math.pi)
       phi1.append(i*10.0/180.0*math.pi)
       for j in range(0,36):
           score_dih.append(0.0)
   # open file

   
   if not mix:

     f = open(filename, 'r')   
     for line in f.readlines():
       riga=(line.strip()).split()
       if (len(riga)==4 and riga[0]==string):
         ii=int(float(riga[1])/10.0)
         jj=int(float(riga[2])/10.0)
         score_dih[ii*len(phi0)+jj]=-KT_CAFF_*math.log(float(riga[3]))
     f.close()        

   if mix:
     #mix random coil and native secondary structure
     counts=[]   
     for i in range(0,36):
       for j in range(0,36):
           counts.append(1.0) 
             
     f = open(filename, 'r')  
     for line in f.readlines():
       riga=(line.strip()).split()
       if (len(riga)==4 and riga[0]==string):
         ii=int(float(riga[1])/10.0)
         jj=int(float(riga[2])/10.0)
         counts[ii*len(phi0)+jj]+=NATIVENESS_*float(riga[3])
       if (len(riga)==4 and riga[0]=="-----"):
         ii=int(float(riga[1])/10.0)
         jj=int(float(riga[2])/10.0)
         counts[ii*len(phi0)+jj]+=(1.0-NATIVENESS_)*float(riga[3])

     f.close()

     for i in range(len(counts)):
        score_dih[i]=-KT_CAFF_*math.log(counts[i])
   

   return [phi0,phi1,score_dih]

def read_potential_angle(filename,string,mix=False):
# read potentials for angles
   score_ang=[]
   psi=[]
   for i in range(0,180):
       psi.append(i/180.0*math.pi)
       score_ang.append(0.0)
   # read file

   
   if not mix:

     f = open(filename, 'r')   
     for line in f.readlines():
       riga=(line.strip()).split()
       if (len(riga)==3 and riga[0]==string):
         ii=int(riga[1])
         score_ang[ii]=-KT_CAFF_*math.log(float(riga[2]))
     f.close()
     
   if mix:
     #mix random coil and native secondary structure
     counts=[]
     for i in range(0,180):
        counts.append(1.0)     

     f = open(filename, 'r')     
     for line in f.readlines():
       riga=(line.strip()).split()
       if (len(riga)==3 and riga[0]==string):
         ii=int(riga[1])
         counts[ii]+=NATIVENESS_*float(riga[2])
       if (len(riga)==3 and riga[0]=="---"): 
         ii=int(riga[1])
         counts[ii]+=(1.0-NATIVENESS_)*float(riga[2])
     f.close()
     

     for i in range(0,180):
        score_ang[i]=-KT_CAFF_*math.log(counts[i])
        
        

   return [psi,score_ang]


def get_CA_force_field(chain, resrange, dihe_dict, ang_dict, do_mix):
    rslist=[]
    pairslist=[]
    # add bonds
    for res in range(resrange[0],resrange[1]): 
        ps=[]
        for delta in range(0,2):
            s=IMP.atom.Selection(chain, residue_index=res+delta, atom_type=IMP.atom.AT_CA)
            ps.append(s.get_selected_particles()[0]) 
        pairslist.append(IMP.ParticlePair(ps[0],ps[1]))
        pairslist.append(IMP.ParticlePair(ps[1],ps[0]))
        br=get_distance_restraint(ps[0],ps[1],3.78,416.0)
        br.set_name('Bond_restraint')
        rslist.append(br)
    # add dihedrals
    for res in range(resrange[0],resrange[1]+1):
        if res not in dihe_dict: continue
        # get the appropriate parameters
        [phi0,phi1,score_dih]=read_potential_dihedral(IMP.isd2.get_data_path("CADihedralRestraint.dat"),dihe_dict[res],do_mix[res])
        # get the particles
        ps=[]
        for delta in range(-2,+3):
            s=IMP.atom.Selection(chain, residue_index=res+delta, atom_type=IMP.atom.AT_CA)
            ps.append(s.get_selected_particles()[0])
        pairslist.append(IMP.ParticlePair(ps[0],ps[3]))
        pairslist.append(IMP.ParticlePair(ps[3],ps[0]))
        pairslist.append(IMP.ParticlePair(ps[1],ps[4]))
        pairslist.append(IMP.ParticlePair(ps[4],ps[1]))
        dr=IMP.isd2.CADihedralRestraint(ps[0],ps[1],ps[2],ps[3],ps[4],phi0,phi1,score_dih)
        dr.set_name('Dihedral restraint')
        rslist.append(dr) 
    # add angles 
    for res in range(resrange[0],resrange[1]+1):
        if res not in ang_dict: continue
        # get the appropriate parameters
        [psi,score_ang]=read_potential_angle(IMP.isd2.get_data_path("CAAngleRestraint.dat"),ang_dict[res],do_mix[res])
        # get the particles
        ps=[]
        for delta in range(-1,+2):
            s=IMP.atom.Selection(chain, residue_index=res+delta, atom_type=IMP.atom.AT_CA)
            ps.append(s.get_selected_particles()[0])
        pairslist.append(IMP.ParticlePair(ps[0],ps[2]))
        pairslist.append(IMP.ParticlePair(ps[2],ps[0]))
        dr=IMP.isd2.CAAngleRestraint(ps[0],ps[1],ps[2],psi,score_ang)
        dr.set_name('Angle restraint')
        rslist.append(dr) 
    return (rslist,pairslist)

def get_homology_restraint(model, maps, sigmaG, cutoff):
    rset=IMP.RestraintSet('Homology_Restraint')
    #dictionaries
    kappas={}
    omegas={}
    dist_dict={}

    for key in maps:
    # define omega and lower and upper bounds
        omega=IMP.isd2.Scale.setup_particle(IMP.Particle(m),1.0)
        omega.set_lower(1.0)
        omega.set_upper(10000.0)
        m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,omega))
        omega.set_is_optimized(omega.get_nuisance_key(),True)
        omegas[key]=omega

    # get all the particles in the model
    ps_model=IMP.atom.get_leaves(model)

    # cycle on all pairs
    for i in range(0,len(ps_model)-1):
        p0=ps_model[i]
        for j in range(i+1,len(ps_model)):
            p1=ps_model[j]
            
            # particles belonging to the same rigid body should not be restrained
            if(IMP.core.RigidMember.particle_is_instance(p0) and 
               IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember(p0).get_rigid_body() == IMP.core.RigidMember(p1).get_rigid_body()): continue

            # cycle over all the dictionaries
            npair=0 
            nclose=0
            for key in maps:
                mk=maps[key]
                if (p0 in mk and p1 in mk): 
                 # get distance in the template 
                 dist0=IMP.core.get_distance(IMP.core.XYZ(mk[p0]),IMP.core.XYZ(mk[p1]))
                 # nclose is the number of templates in which the pair distance is lower than cutoff
                 if(dist0<cutoff): nclose+=1
                 # npair is the number of templates in which the pair exists
                 npair+=1
            # if the pair exists in all the templates and it is lower than the cutoff at least in one, add restraint
            if(npair==len(maps) and nclose>0):
              for key in maps:
                mk=maps[key] 
                dist_dict[(key,p0,p1)]=IMP.core.get_distance(IMP.core.XYZ(mk[p0]),IMP.core.XYZ(mk[p1]))
                
              # define kappa and lower and upper bounds
              kappa=IMP.isd2.Scale.setup_particle(IMP.Particle(m),1.0)
              kappa.set_lower(0.0)
              kappa.set_upper(float(npair-1))
              m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,kappa))
              if(npair>1): kappa.set_is_optimized(kappa.get_nuisance_key(),True)
              else: kappa.set_is_optimized(kappa.get_nuisance_key(),False)
              r0=IMP.atom.Residue(IMP.atom.Atom(p0).get_parent()).get_index()
              r1=IMP.atom.Residue(IMP.atom.Atom(p1).get_parent()).get_index() 
              kappas[str(r0)+"-"+str(r1)]=kappa
              lnar=IMP.isd2.LogNormalAmbiguousRestraint(p0,p1,kappa,sigmaG)
              for key in maps:
                   lnar.add_contribution(dist_dict[(key,p0,p1)],omegas[key])
              rset.add_restraint(lnar)

    return kappas,omegas,rset,dist_dict

def get_elastic_restraint(model, residlist, sigmaG, cutoff):
    rslist=[]

    # define omega and lower and upper bounds
    omega=IMP.isd2.Scale.setup_particle(IMP.Particle(m),1.0)
    omega.set_lower(1.0)
    omega.set_upper(10000.0)
    m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,omega))
    omega.set_is_optimized(omega.get_nuisance_key(),True)

    # get all the particles in the model
    s=IMP.atom.Selection(model, residue_indexes=residlist)
    ps_model=s.get_selected_particles()

    # cycle on all pairs
    for i in range(0,len(ps_model)-1):
        p0=ps_model[i]
        for j in range(i+1,len(ps_model)):
            p1=ps_model[j]
            
            # particles belonging to the same rigid body should not be restrained
            if(IMP.core.RigidMember.particle_is_instance(p0) and 
               IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember(p0).get_rigid_body() == IMP.core.RigidMember(p1).get_rigid_body()): continue

            # get distance in the template 
            dist0=IMP.core.get_distance(IMP.core.XYZ(p0),IMP.core.XYZ(p1))
            if(dist0<cutoff):
              # define kappa and lower and upper bounds
              kappa=IMP.isd2.Scale.setup_particle(IMP.Particle(m),0.0)
              kappa.set_lower(0.0)
              kappa.set_upper(0.0)
              m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,kappa))
              kappa.set_is_optimized(kappa.get_nuisance_key(),False)
              lnar=IMP.isd2.LogNormalAmbiguousRestraint(p0,p1,kappa,sigmaG)
              lnar.add_contribution(dist0,omega)
              rslist.append(lnar)

    return omega,rslist

def get_kink_restraint(rbpairs, kappa):
    rset=IMP.RestraintSet('Kink_Restraint')
    for rbs in rbpairs:
        h=IMP.core.HarmonicWell([10.0/180.0*math.pi,40.0/180.0*math.pi],kappa)
        kps=IMP.membrane.KinkPairScore(h)
        pr=IMP.core.PairRestraint(kps,IMP.ParticlePair(rbs[0],rbs[1]))
        rset.add_restraint(pr)
    return rset
    
def get_prior(sigmas):
    rslist=[]
    for sigma in sigmas:
        rslist.append(IMP.isd2.JeffreysRestraint(sigmas[sigma]))
    return rslist

def get_rb_movers(rblist,rmax,tmax):
    mvs=[]
    for rb in rblist:
        mv= IMP.core.RigidBodyMover(rblist[rb], rmax, tmax)
        mvs.append(mv)
    return mvs

def get_ball_movers(ps,tmax):
    mvs=[]
    for p in ps: 
        mv= IMP.core.BallMover([p], tmax)
        mvs.append(mv)
    return mvs

def get_nuisance_movers(omegas,maxstep):
    mvs=[]
    for omega in omegas:
        mvs.append(IMP.core.NormalMover([omegas[omega]],IMP.FloatKeys([IMP.FloatKey("nuisance")]),maxstep))
    return mvs

def get_drms(align,dist_dict):
    drms={}
    npairs={} 
    for key in align:
        drms[key]=0.0    
        npairs[key]=0.0

    for key in dist_dict:
        k0=key[0]
        p0=key[1]
        p1=key[2]
        dist=IMP.core.get_distance(IMP.core.XYZ(p0),IMP.core.XYZ(p1)) 
        dist0=dist_dict[key]
        drms[key[0]]+=(dist-dist0)*(dist-dist0)
        npairs[key[0]]+=1.0
   
    for key in drms:
        drms[key]=math.sqrt(drms[key]/npairs[key])

    return drms

def simulated_annealing(nuisances,signs,istep,ncold,nhot):
     if istep%(ncold+nhot)< ncold:
        value=0.0
     else:
        a=float(istep%(ncold+nhot)-ncold)
        value=1.0/4.0*((1.0-math.cos(2*math.pi*a/float(nhot))))**2
# nuisances is a list of nuisance parameters,
# signs is a list of sign for the direction of the annealing
     for i,nus in enumerate(nuisances):
         nus.set_scale(0.5*(1.0+signs[i])*nus.get_lower()+0.5*(1.0-signs[i])*nus.get_upper()+signs[i]*(nus.get_upper()-nus.get_lower())*value)

def get_grid(gmin,gmax,ngrid,boundaries):                                        
    grid=[]                                                                      
    dx = ( gmax - gmin ) / float(ngrid)                                          
    for i in range(0,ngrid+1):                                                   
        if(not boundaries and i==0): continue                                    
        if(not boundaries and i==ngrid): continue                                
        grid.append( gmin + float(i) * dx )                                      
    return grid                                                                  
                                                                                 
def get_log_grid(gmin,gmax,ngrid):                                               
   grid=[]                                                                       
   for i in range(0,ngrid+1):                                                    
       grid.append( gmin*math.exp(float(i)/ngrid*math.log(gmax/gmin)) )                    
   return grid

def get_crosslink_restraint(protein_copies, filename):
    rset=IMP.RestraintSet('Cysteine_Crosslink')
    # dictionaries
    expdict={}
    epsilons={}

    # beta
    beta=IMP.isd2.Scale.setup_particle(IMP.Particle(m),BETA_)
    beta.set_lower(BETA_)
    beta.set_upper(BETA_)
    m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,beta))
    beta.set_is_optimized(beta.get_nuisance_key(),False)

    # sigma 
    sigma=IMP.isd2.Scale.setup_particle(IMP.Particle(m),1.0)
    sigma.set_lower(0.0)
    sigma.set_upper(1.0)
    m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,sigma))
    sigma.set_is_optimized(sigma.get_nuisance_key(),False)

    # upper boundary dictionary
    upperb={}
    upperb["16-41"]=1.0-0.628809
    upperb["42-65"]=1.0-0.922462
    upperb["185-226"]=1.0-0.693342
 
    # epsilon particles 
    for id in ["16-41","42-65","185-226"]:
        epsilonscale=IMP.isd2.Scale.setup_particle(IMP.Particle(m),0.01)
        epsilonscale.set_lower(0.01)
        epsilonscale.set_upper(upperb[id])
        m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,epsilonscale))
        epsilonscale.set_is_optimized(epsilonscale.get_nuisance_key(),True)
        epsilons[id]=epsilonscale
   
    # population particle 
    pw=IMP.Particle(m)                                                               
    weight=IMP.isd2.Weight.setup_particle(pw)                                        
    weight.set_weights_are_optimized(True)

    # create grids needed by CrossLinkData
    dist_grid=get_grid(0.0,25.0, NGRID_, False)  
    omega1_grid=get_log_grid(1.0, 1000.0, 50)
    sigma_grid=[1.0]
    # read PMF from file
    xlpot=open("cysteine_FES.dat")                                                   
    pot_x_grid=[]                                                                      
    pot_value_grid=[]                                                                   
    for line in xlpot:                                                               
        t=line.split()                                                               
        if t[0][0]!="#":                                                             
           x = float(t[0])*10.0
           pot_x_grid.append(x)
           pot_value_grid.append(float(t[1])/4.0/math.pi/x/x)
    # CrossLinkMSData 
    crossdata=IMP.isd2.CrossLinkData(dist_grid,omega1_grid,sigma_grid,pot_x_grid,pot_value_grid,10.0,20.0)

    # create grids needed by CysteineCrossLinkData
    fmod_grid=get_grid(0.0, 1.0, 300, True)
    omega2_grid=get_log_grid(0.001, 10000.0, 100) 

    # open file with cross-link data
    f=open(filename,"r")
    # read file
    for line in f.readlines():
        riga=(line.strip()).split()
        resid1=int(riga[0])
        chain1=riga[1]
        resid2=int(riga[2])
        chain2=riga[3]
        fexp=float(riga[4])

        # CysteineCrossLinkData
        ccldata=IMP.isd2.CysteineCrossLinkData(fexp,fmod_grid,omega2_grid,[BETA_])

        # get upperb id
        id=""
        if( 16 <= resid1 <= 41):   id="16-41"
        if( 42 <= resid1 <= 65):   id="42-65"
        if(185 <= resid1 <=226): id="185-226"

        ccl=IMP.isd2.CysteineCrossLinkRestraint(beta,sigma,epsilons[id],pw,crossdata,ccldata)

        for i,prot in enumerate(protein_copies):
            s1=IMP.atom.Selection(prot, chains=chain1, residue_index=resid1, atom_type=IMP.atom.AT_CA)
            p1=s1.get_selected_particles()[0]

            s2=IMP.atom.Selection(prot, chains=chain2, residue_index=resid2, atom_type=IMP.atom.AT_CA)
            p2=s2.get_selected_particles()[0]

            ccl.add_contribution(p1,p2)

        ccl.set_name(str(resid1))
        rset.add_restraint(ccl) 

        # build dictionaries
        expdict[str(resid1)]=fexp

    return rset,expdict,beta,epsilons,pw

def get_Ez_potential(protein,boundaries): 
    rset=IMP.RestraintSet('Ez_potential')
    ps=[]
    for boundary in boundaries: 
        s=IMP.atom.Selection(protein, residue_indexes=boundary)
        ps+=s.get_selected_particles()
    ez=IMP.membrane.EzRestraint(ps)
    rset.add_restraint(ez)
    return rset 

def initialize_coordinates(protein):
    rb=IMP.atom.create_rigid_body(protein)
    rbcoord=rb.get_coordinates()
    rot=IMP.algebra.get_identity_rotation_3d()
    tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
    trans=tmptrans.get_inverse()
    IMP.core.transform(rb,trans)
    rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(IMP.algebra.Transformation3D
       (IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0,1,0),-math.pi/2.0),
        IMP.algebra.Vector3D(0,0,+14.0))))
    IMP.core.RigidBody.teardown_particle(rb)
    m.remove_particle(rb)
 
def get_layer_restraint(protein,resid,zrange,kappa):
    rset=IMP.RestraintSet('Layer_restraint')
    lsc=IMP.container.ListSingletonContainer(m)
    s=IMP.atom.Selection(protein, residue_index=resid, atom_type=IMP.atom.AT_CA)
    lsc.add_particles(s.get_selected_particles())
    hw=IMP.core.HarmonicWell(zrange,kappa) 
    asc=IMP.core.AttributeSingletonScore(hw,IMP.FloatKey("z"))
    sr=IMP.container.SingletonsRestraint(asc, lsc)
    rset.add_restraint(sr)
    return rset

# read initial model
PHOQ=[]
#filelist=["model_A.pdb","model_B.pdb"]
filelist=["PHOQ_PERI_initial_model.pdb","PHOQ_PERI_initial_model.pdb","PHOQ_PERI_initial_model.pdb"]
for i in range(0,NCOPIES_):
    PHOQ.append(load_pdb("PHOQ copy "+str(i),filelist[i]))

# initialize system position
for i in range(0,NCOPIES_): initialize_coordinates(PHOQ[i])

# create rigid bodies
rblist={}
fixres=[(13,41),(45,184),(194,205),(208,217),(220,233),(245,265)]
for i in range(0,NCOPIES_):
    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
        for res in fixres:
            rb=create_rigid_body(chain,res[0],res[1])
            rblist[IMP.atom.Chain(chain).get_id()+"_"+str(res[0])+"-"+str(res[1])+"::"+str(i)]=rb

# Create GLOBAL restraints set
global_rset={}

# Excluded Volume
for i in range(0,NCOPIES_):
    global_rset['Excluded_Volume::'+str(i)]=get_excluded_volume(PHOQ[i],1.0)

# Link rigid and floppy parts
floppyres=[(41,45),(184,194),(205,208),(217,220),(233,245)]

# CA force field dictionaries
dihe_dict={}; ang_dict={}; do_mix={}
# 41-45: from 42 to 44 coiled-coil, helical otherwise
dihe_dict[41]="HHH--";     ang_dict[41]="HH-";   do_mix[41]=False
dihe_dict[42]="HH---";     ang_dict[42]="H--";   do_mix[42]=False
dihe_dict[43]="H---H";     ang_dict[43]="---";   do_mix[43]=False
dihe_dict[44]="---HH";     ang_dict[44]="--H";   do_mix[44]=False
dihe_dict[45]="--HHH";     ang_dict[45]="-HH";   do_mix[45]=False
# 184-194: from 186 to 189 coiled-coil, helical otherwise 
dihe_dict[184]="HHHH-";    ang_dict[184]="HHH";  do_mix[184]=False
dihe_dict[185]="HHH--";    ang_dict[185]="HH-";  do_mix[185]=False
dihe_dict[186]="HH---";    ang_dict[186]="H--";  do_mix[186]=False
dihe_dict[187]="H----";    ang_dict[187]="---";  do_mix[187]=False
dihe_dict[188]="----H";    ang_dict[188]="---";  do_mix[188]=False
dihe_dict[189]="---HH";    ang_dict[189]="--H";  do_mix[189]=False
dihe_dict[190]="--HHH";    ang_dict[190]="-HH";  do_mix[190]=False
dihe_dict[191]="-HHHH";    ang_dict[191]="HHH";  do_mix[191]=False
dihe_dict[192]="HHHHH";    ang_dict[192]="HHH";  do_mix[192]=False
dihe_dict[193]="HHHHH";    ang_dict[193]="HHH";  do_mix[193]=False
dihe_dict[194]="HHHHH";    ang_dict[194]="HHH";  do_mix[194]=False
# 205-208: PRO 208 kink 
dihe_dict[205]="HHHHH";    ang_dict[205]="HHH"; do_mix[205]=True
dihe_dict[206]="HHHHH";    ang_dict[206]="HHH"; do_mix[206]=True
dihe_dict[207]="HHHHH";    ang_dict[207]="HHH"; do_mix[207]=True
dihe_dict[208]="HHHHH";    ang_dict[208]="HHH"; do_mix[208]=True
# 217-220: PRO 220 kink
dihe_dict[217]="HHHHH";    ang_dict[217]="HHH"; do_mix[217]=True
dihe_dict[218]="HHHHH";    ang_dict[218]="HHH"; do_mix[218]=True
dihe_dict[219]="HHHHH";    ang_dict[219]="HHH"; do_mix[219]=True
dihe_dict[220]="HHHHH";    ang_dict[220]="HHH"; do_mix[220]=True
# random coil in the hamp domain
for res in range(233,246):
    dihe_dict[res]="-----";    ang_dict[res]="---";  do_mix[res]=False

# dictionary of ListPairContainers for bonded interactions
excvol_filter={}

for i in range(0,NCOPIES_):
    global_rset['CAForceField::'+str(i)]=IMP.RestraintSet()
    excvol_filter['Excluded_Volume::'+str(i)]=IMP.container.ListPairContainer(m)
    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
        for res in floppyres:
            (rstlist,pairslist)=get_CA_force_field(chain,res,dihe_dict,ang_dict,do_mix)
            excvol_filter['Excluded_Volume::'+str(i)].add_particle_pairs(pairslist)
            global_rset['CAForceField::'+str(i)].add_restraints(rstlist)

# create excluded volume filter for bonded pairs
for key in excvol_filter:
    # get the ListPairContainer from dictionary
    lpc=excvol_filter[key]
    # create filter
    icpf=IMP.container.InContainerPairFilter(lpc)
    # get restraint from global dictionary
    rst=global_rset[key].get_restraints()[0]
    # add filter
    IMP.core.ExcludedVolumeRestraint.get_from(rst).add_pair_filter(icpf)

# list of ISD particles (needed to write to RMF)
ISD_particles=[]

# sigmaG is the same for all homology and elastic restraints
sigmaG=IMP.isd2.Scale.setup_particle(IMP.Particle(m),SIGMAG_)
sigmaG.set_lower(SIGMAG_)
sigmaG.set_upper(SIGMAG_)
m.add_score_state(IMP.core.SingletonConstraint(IMP.isd2.NuisanceRangeModifier(),None,sigmaG))
sigmaG.set_is_optimized(sigmaG.get_nuisance_key(),False)
ISD_particles.append(sigmaG) # add to list of ISD particles for rmf I/O

# Elastic Network
#elastic_omegas_copies=[]
#elasticres=[(37,49),(181,198),(213,225)]
#elasticres=[(37,49),(181,198),(203,213),(213,225)]
#elasticres=[(37,49),(203,213),(213,225)]

#for i in range(0,NCOPIES_):
#    elastic_omegas={}
#    global_rset['Elastic_Network::'+str(i)]=IMP.RestraintSet()
#    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
#        for res in elasticres:
#            (omega,rs)=get_elastic_restraint(chain, [res], sigmaG, ELASTIC_CUTOFF_)
#            global_rset['Elastic_Network::'+str(i)].add_restraints(rs)
#            elastic_omegas[IMP.atom.Chain(chain).get_id()+"_"+str(res[0])+"-"+str(res[1])]=omega
#            ISD_particles.append(omega) # add to list of ISD particles for rmf I/O
#    elastic_omegas_copies.append(elastic_omegas)

# Homology restraints
# list for multiple copies sampling
hamp_align_copies=[]
hamp_kappas_copies=[]
hamp_omegas_copies=[]
hamp_dist_dict_copies=[]

# HAMP DOMAIN
# 1) 3ZRX
#template1=load_pdb("template1","3ZRX.pdb")
# 2) 2Y20_A_B 
template2=load_pdb("template2","2Y20_A_B.pdb")

# PERIPLASMIC
#peri_align_copies=[]
#peri_kappas_copies=[]
#peri_omegas_copies=[]
#peri_dist_dict_copies=[]
#template3=load_pdb("template3","3BQ8.pdb")

for i in range(0,NCOPIES_):
    hamp_align={}
    #hamp_align['3ZRX']=read_alignments("PHOQ-3ZRX.align",PHOQ[i],template1)
    hamp_align['2Y20']=read_alignments("PHOQ-2Y20.align",PHOQ[i],template2)
    hamp_align_copies.append(hamp_align)

    (hamp_kappas,hamp_omegas,rset,hamp_dist_dict)=get_homology_restraint(PHOQ[i],hamp_align,sigmaG,HOMOLOGY_CUTOFF_)
    hamp_kappas_copies.append(hamp_kappas)
    hamp_omegas_copies.append(hamp_omegas)
    hamp_dist_dict_copies.append(hamp_dist_dict)
     
    global_rset['Homology_Hamp::'+str(i)]=rset
    # add to list of ISD particles for rmf I/O
    for key in hamp_kappas:
        ISD_particles.append(hamp_kappas[key])
    for key in hamp_omegas:
        ISD_particles.append(hamp_omegas[key])   

    #peri_align={}
    #peri_align['3BQ8']=read_alignments("PHOQ-3BQ8.align",PHOQ[i],template3)
    #peri_align_copies.append(peri_align)

    #(peri_kappas,peri_omegas,rset,peri_dist_dict)=get_homology_restraint(PHOQ[i],peri_align,sigmaG,HOMOLOGY_CUTOFF_)
    #peri_kappas_copies.append(peri_kappas)
    #peri_omegas_copies.append(peri_omegas)
    #peri_dist_dict_copies.append(peri_dist_dict)

    #global_rset['Homology_Perii::'+str(i)]=rset
    # add to list of ISD particles for rmf I/O
    #for key in peri_kappas:
    #    ISD_particles.append(peri_kappas[key])
    #for key in peri_omegas:
    #    ISD_particles.append(peri_omegas[key]) 


# Cross-link restraint
(global_rset['Cross-link'],crosslink_expdict,crosslink_beta,crosslink_epsilons,crosslink_pw)=get_crosslink_restraint(PHOQ,"expcrosslink_95_percent.dat")
# add to list of ISD particles for rmf I/O
ISD_particles.append(crosslink_beta)
ISD_particles.append(crosslink_pw)
for key in crosslink_epsilons:
    ISD_particles.append(crosslink_epsilons[key])

# Jeffrey Prior for omegas
for i in range(0,NCOPIES_):
   global_rset['Prior::'+str(i)]=IMP.RestraintSet()
   global_rset['Prior::'+str(i)].add_restraints(get_prior(hamp_omegas_copies[i]))
   #global_rset['Prior::'+str(i)].add_restraints(get_prior(peri_omegas_copies[i]))
   #global_rset['Prior::'+str(i)].add_restraints(get_prior(elastic_omegas_copies[i]))

# Boundaries for Weight
global_rset['Weight Boundaries']=IMP.isd2.WeightRestraint(crosslink_pw,0.2,0.8,100000.)

# Ez potential
#TM_regions=[[(13,45)],[(194,221)]]
TM_regions=[[(17,37)],[(195,215)]]
for i in range(0,NCOPIES_):
    global_rset['Ez_potential::'+str(i)]=get_Ez_potential(PHOQ[i],TM_regions)

# Layer restraint
N_terminus=17
layer=(-17.0,-13.0)
for i in range(0,NCOPIES_):
    global_rset['Layer_restraint::'+str(i)]=get_layer_restraint(PHOQ[i],N_terminus,layer,100.0*KAPPA_)

# kink restraint
#for i in range(0,NCOPIES_):
#    hpairs=[]
#    hpairs.append([rblist["A_194-205::"+str(i)],rblist["A_208-217::"+str(i)]])
#    hpairs.append([rblist["B_194-205::"+str(i)],rblist["B_208-217::"+str(i)]])
#    global_rset['Kink::'+str(i)]=get_kink_restraint(hpairs, 100.0*KAPPA_)

# add restraints to model
for rset in global_rset:
    m.add_restraint(global_rset[rset])

# Set coordinates as not optimized
for i in range(0,NCOPIES_):
    for atom in IMP.atom.get_leaves(PHOQ[i]): 
        IMP.core.XYZR(atom).set_coordinates_are_optimized(True)

# Movers
mvs=[]
# Rigid Body Movers
mvs+=get_rb_movers(rblist,MAXTRANS_,MAXROT_)
# Particle Movers
ps=[]
for i in range(0,NCOPIES_):
    for p in IMP.atom.get_leaves(PHOQ[i]):
        if(IMP.core.RigidMember.particle_is_instance(p)==False):
          ps.append(p)
mvs+=get_ball_movers(ps,MAXTRANS_)

# Movers for omegas 
for i in range(0,NCOPIES_):
    mvs+=get_nuisance_movers(hamp_omegas_copies[i],MAXOMEGA_)
    #mvs+=get_nuisance_movers(peri_omegas_copies[i],MAXOMEGA_)
    #mvs+=get_nuisance_movers(elastic_omegas_copies[i],MAXOMEGA_)

# Movers for epsilons 
mvs+=get_nuisance_movers(crosslink_epsilons,MAXEPSILON_) 

# Movers for weight
mvs.append(IMP.isd2.WeightMover(crosslink_pw,MAXWEIGHT_))

# SerialMover
smv=IMP.core.SerialMover(mvs)

# Monte Carlo with Wte
#mc=IMP.membrane.MonteCarloWithWte(m,EMIN_,EMAX_,SIGMA_,GAMMA_,W0_)
mc=IMP.core.MonteCarlo(m)
mc.set_return_best(False)
mc.set_kt(tempval)
mc.add_mover(smv)

# add bias parameters
#mybias=mc.get_bias_asfloats()
#rem.set_my_parameter("bias",mybias)

# open rmf for writing
rh = RMF.create_rmf_file("traj"+str(myindex)+".rmf")
# add coordinates
for i in range(0,NCOPIES_): IMP.rmf.add_hierarchy(rh, PHOQ[i])
# add ISD particles
IMP.rmf.add_particles(rh, ISD_particles)
# add other information
my_kc=rh.add_category("my data")
score_key=rh.add_float_key(my_kc,"my score",True)
cross_key=rh.add_float_key(my_kc,"my cross-link score",True)
#bias_key=rh.add_float_key(my_kc,"my bias",True)
index_key=rh.add_int_key(my_kc,"my index",True)

# open log file
log=open('log'+str(myindex),'w')

# Sampling (or reloading)
for istep in range(0,NITER_):

    # Monte Carlo
    mc.optimize(NOPT_)
 
    #for i in range(0,len(mvs)):
    #    print i,smv.get_acceptance(i)

    # time for printout
    if(istep%W_STRIDE==0):

        hamp_drms_copies=[]
        #peri_drms_copies=[]
        for i in range(0,NCOPIES_):
        # get DRMS for hamp and peri domain
            hamp_drms_copies.append(get_drms(hamp_align_copies[i],hamp_dist_dict_copies[i]))
            #peri_drms_copies.append(get_drms(peri_align_copies[i],peri_dist_dict_copies[i]))
 
        # get total score, bias potential and cross-link score for RMF
        myscore=m.evaluate(False)
        #mybias=mc.get_bias(myscore)
        mycross=global_rset['Cross-link'].evaluate(False)

        # get weights                                                            
        ww=IMP.isd2.Weight(crosslink_pw).get_weights()
 
        # prepare printout
        s0="%12d " % (istep)
        s00="%12d " % (myindex)
        #s000="%12.6f " % (mybias)   
        s0000="%12.6f " % (mc.get_kt())
        s1="%12.6f " % (mc.get_number_of_forward_steps()/float(NOPT_))
        s2=' '.join(["%5s %12.6f " % (kkey,global_rset[kkey].evaluate(False)) for kkey in global_rset])
        #s6=' '.join(["%5s %12.6f " % (kkey,sigmas[kkey].get_scale()) for kkey in sigmas])
        #s7=' '.join(["%5s %12.6f " % (kkey,betas[kkey].get_scale()) for kkey in betas])
        s9=' '.join(["%5s %12.6f " % (IMP.isd2.CysteineCrossLinkRestraint.get_from(rst).get_name(),IMP.isd2.CysteineCrossLinkRestraint.get_from(rst).get_model_frequency()) for rst in global_rset['Cross-link'].get_restraints()])
        s10=' '.join(["%5s %12.6f " % (kkey,crosslink_expdict[kkey]) for kkey in crosslink_expdict])
        s11=' '.join(["%5s %12.6f " % (IMP.isd2.CysteineCrossLinkRestraint.get_from(rst).get_name(),IMP.isd2.CysteineCrossLinkRestraint.get_from(rst).get_standard_error()) for rst in global_rset['Cross-link'].get_restraints()]) 
        s14="Total_energy %lf" % (myscore)
        s20=' '.join(["%5s %12.6f " % (kkey,crosslink_epsilons[kkey].get_scale()) for kkey in crosslink_epsilons])
        s30=' '.join(["%5d %12.6f " % (i,ww[i]) for i in range(0,ww.get_dimension())]) 

        # printout copy-dependent
        s3=[]; s4=[]; s5=[]; s12=[]; s13=[];
        for i in range(0,NCOPIES_):
            s3.append(' '.join(["hamp %5s %12.6f " % (kkey,hamp_omegas_copies[i][kkey].get_scale()) for kkey in hamp_omegas_copies[i]]))
            #s4.append(' '.join(["peri %5s %12.6f " % (kkey,peri_omegas_copies[i][kkey].get_scale()) for kkey in peri_omegas_copies[i]]))
            #s5.append(' '.join(["%5s %12.6f " % (kkey,elastic_omegas_copies[i][kkey].get_scale()) for kkey in elastic_omegas_copies[i]]))
            s12.append(' '.join(["hamp %5s %12.6f " % (kkey,hamp_drms_copies[i][kkey]) for kkey in hamp_drms_copies[i]]))
            #s13.append(' '.join(["peri %5s %12.6f " % (kkey,peri_drms_copies[i][kkey]) for kkey in peri_drms_copies[i]]))
  
 
        log.write("REM_Index                 %s" % s0+s00             )
        log.write("\n")
        log.write("Acceptance                %s" % s0+" "+s1          )
        log.write("\n")
        log.write("Scores                    %s" % s0+" "+s14+" "+s2  )
        log.write("\n")
        log.write("Temperature               %s" % s0+" "+s0000       )
        log.write("\n")
        #log.write("Bias_potential            %s" % s0+" "+s000        )
        #log.write("\n")
        for i in range(0,NCOPIES_):
            log.write("Homology_Omegas           %s" % s0+" copyID_"+str(i)+" "+s3[i]) #+" "+s4[i]   )
            log.write("\n")
        #for i in range(0,NCOPIES_):
        #    log.write("Elastic_Omegas            %s" % s0+" copyID_"+str(i)+" "+s5[i]          )
        #    log.write("\n")

        log.write("Crosslink_nu              %s" % s0+" "+s9          )
        log.write("\n")
        log.write("Crosslink_exp_data        %s" % s0+" "+s10         )
        log.write("\n")
        log.write("Crosslink_standard_error  %s" % s0+" "+s11         )
        log.write("\n") 
        #log.write("Sigmas                    %s" % s0+" "+s6          )
        #log.write("\n")
        #log.write("Betas                     %s" % s0+" "+s7          )
        #log.write("\n")
        log.write("Crosslink_epsilons        %s" % s0+" "+s20         )
        log.write("\n")
        log.write("Crosslink_weights         %s" % s0+" "+s30         )
        log.write("\n")

        for i in range(0,NCOPIES_):
            log.write("Drms                      %s" % s0+" copyID_"+str(i)+" "+s12[i]) #+" "+s13[i] )
            log.write("\n")
 
        # print all information of the current frame to rmf
        IMP.rmf.save_frame(rh,istep/W_STRIDE)
        # and other useful data
        (rh.get_root_node()).set_value(score_key,myscore,istep/W_STRIDE)
        (rh.get_root_node()).set_value(cross_key,mycross,istep/W_STRIDE)
        #(rh.get_root_node()).set_value(bias_key,mybias,istep/W_STRIDE)
        (rh.get_root_node()).set_value(index_key,myindex,istep/W_STRIDE)


    # time for an exchange
    score=m.evaluate(False)
    # get my index
    myindex=rem.get_my_index()
    # get my parameters 
    mytemp=rem.get_my_parameter("temp")[0]
    #mybias=mc.get_bias_asfloats()
    #rem.set_my_parameter("bias",mybias)
    # calculate model score
    myscore=score / mytemp

    # get friend's index 
    findex=rem.get_friend_index(istep)
    # get friend's parameters 
    ftemp=rem.get_friend_parameter("temp",findex)[0]
    # calculate model score
    fscore=score / ftemp

    # try exchange
    flag=rem.do_exchange(myscore,fscore,findex)
    # if accepted, change temperature 
    if (flag==True):
       mc.set_kt(ftemp)

# close all files
log.close()
rh.flush()
rh=RMF.FileHandle()
