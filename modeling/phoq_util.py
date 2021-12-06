from __future__ import print_function
import IMP.atom
import math

try:
    from IMP.mpi import ReplicaExchange
except ImportError:
    from IMP.pmi.samplers import _SerialReplicaExchange as ReplicaExchange

# MIXING PARAMETER FOR CA FF
NATIVENESS_=0.3
KT_CAFF_=1.0

# CROSSLINK PARA
# typical beta derived from average standard deviation
# of experimental error
BETA_=0.025
SIGMAG_=0.02

NGRID_=1000

def read_alignments(filename,model,temp):
    f=open(filename,"r")
    align={}
    for line in f.readlines():
        riga=(line.strip()).split()

        try:
            s1=IMP.atom.Selection(model, chains=riga[1], residue_index=int(riga[2]), residue_type=IMP.atom.get_residue_type(riga[0]), atom_type=IMP.atom.AT_CA)
            p1=s1.get_selected_particles()[0]
        except:
            print("residue %d of chain %s not found in model" % (int(riga[2]),riga[1]))
            continue
        try:
            s2=IMP.atom.Selection(temp, chains=riga[4], residue_index=int(riga[5]), residue_type=IMP.atom.get_residue_type(riga[3]), atom_type=IMP.atom.AT_CA)
            p2=s2.get_selected_particles()[0]
        except:
            print("residue %d of chain %s not found in template" % (int(riga[5]), riga[4]))
            continue

        align[p1]=p2

    f.close()
    return align

def load_pdb(m, name, filename):
    prot=IMP.atom.read_pdb(filename,m,IMP.atom.CAlphaPDBSelector())
    prot.set_name(name)
    return prot

def create_rigid_body(m, chain,res0,res1):
    s=IMP.atom.Selection(chain, residue_indexes=range(res0,res1+1))
    atoms=[]
    for p in s.get_selected_particles():
        atoms.append(IMP.core.XYZR(p))
    prb=IMP.Particle(m)
    rb=IMP.core.RigidBody.setup_particle(prb,atoms)
    return rb

def get_excluded_volume(m, prot,kappa,excluded_residues=[]):
    rset=IMP.RestraintSet(m, 'Excluded_Volume')
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
    lsc.add(atoms_to_use)
    evr=IMP.core.ExcludedVolumeRestraint(lsc,kappa)
    evr.set_name('Excluded_Volume')
    rset.add_restraint(evr)
    return rset

def get_distance_restraint(m, p0, p1, d0, kappa):
    h=IMP.core.Harmonic(d0,kappa)
    dps=IMP.core.DistancePairScore(h)
    pr=IMP.core.PairRestraint(m, dps,IMP.ParticlePair(p0,p1))
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


def get_CA_force_field(m, chain, resrange, dihe_dict, ang_dict, do_mix):
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
        br=get_distance_restraint(m, ps[0],ps[1],3.78,416.0)
        br.set_name('Bond_restraint')
        rslist.append(br)
    # add dihedrals
    for res in range(resrange[0],resrange[1]+1):
        if res not in dihe_dict: continue
        # get the appropriate parameters
        [phi0,phi1,score_dih]=read_potential_dihedral("../data/CADihedralRestraint.dat",dihe_dict[res],do_mix[res])
        # get the particles
        ps=[]
        for delta in range(-2,+3):
            s=IMP.atom.Selection(chain, residue_index=res+delta, atom_type=IMP.atom.AT_CA)
            ps.append(s.get_selected_particles()[0])
        pairslist.append(IMP.ParticlePair(ps[0],ps[3]))
        pairslist.append(IMP.ParticlePair(ps[3],ps[0]))
        pairslist.append(IMP.ParticlePair(ps[1],ps[4]))
        pairslist.append(IMP.ParticlePair(ps[4],ps[1]))
        dr=IMP.atom.CADihedralRestraint(m, ps[0],ps[1],ps[2],ps[3],ps[4],phi0,phi1,score_dih)
        dr.set_name('Dihedral restraint')
        rslist.append(dr)
    # add angles
    for res in range(resrange[0],resrange[1]+1):
        if res not in ang_dict: continue
        # get the appropriate parameters
        [psi,score_ang]=read_potential_angle("../data/CAAngleRestraint.dat",ang_dict[res],do_mix[res])
        # get the particles
        ps=[]
        for delta in range(-1,+2):
            s=IMP.atom.Selection(chain, residue_index=res+delta, atom_type=IMP.atom.AT_CA)
            ps.append(s.get_selected_particles()[0])
        pairslist.append(IMP.ParticlePair(ps[0],ps[2]))
        pairslist.append(IMP.ParticlePair(ps[2],ps[0]))
        dr=IMP.atom.CAAngleRestraint(m, ps[0],ps[1],ps[2],psi,score_ang)
        dr.set_name('Angle restraint')
        rslist.append(dr)
    return (rslist,pairslist)

def get_homology_restraint(m, model, maps, sigmaG, cutoff):
    rset=IMP.RestraintSet(m, 'Homology_Restraint')
    #dictionaries
    kappas={}
    omegas={}
    dist_dict={}

    for key in maps:
    # define omega and lower and upper bounds
        omega=IMP.isd.Scale.setup_particle(IMP.Particle(m),1.0)
        omega.set_lower(1.0)
        omega.set_upper(10000.0)
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
            if(IMP.core.RigidMember.get_is_setup(p0) and
               IMP.core.RigidMember.get_is_setup(p1) and
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
                kappa=IMP.isd.Scale.setup_particle(IMP.Particle(m),1.0)
                kappa.set_lower(0.0)
                kappa.set_upper(float(npair-1))
                if(npair>1): kappa.set_is_optimized(kappa.get_nuisance_key(),True)
                else: kappa.set_is_optimized(kappa.get_nuisance_key(),False)
                r0=IMP.atom.Residue(IMP.atom.Atom(p0).get_parent()).get_index()
                r1=IMP.atom.Residue(IMP.atom.Atom(p1).get_parent()).get_index()
                kappas[str(r0)+"-"+str(r1)]=kappa
                lnar=IMP.isd.LognormalAmbiguousRestraint(p0,p1,kappa,sigmaG)
                for key in maps:
                    lnar.add_contribution(dist_dict[(key,p0,p1)],omegas[key])
                rset.add_restraint(lnar)

    return kappas,omegas,rset,dist_dict

def get_elastic_restraint(m, model, residlist, sigmaG, cutoff):
    rslist=[]

    # define omega and lower and upper bounds
    omega=IMP.isd.Scale.setup_particle(IMP.Particle(m),1.0)
    omega.set_lower(1.0)
    omega.set_upper(10000.0)
    omega.set_is_optimized(omega.get_nuisance_key(),True)

    # get all the particles in the model
    s=IMP.atom.Selection(model, residue_indexes=range(residlist))
    ps_model=s.get_selected_particles()

    # cycle on all pairs
    for i in range(0,len(ps_model)-1):
        p0=ps_model[i]
        for j in range(i+1,len(ps_model)):
            p1=ps_model[j]

            # particles belonging to the same rigid body should not be restrained
            if(IMP.core.RigidMember.get_is_setup(p0) and
               IMP.core.RigidMember.get_is_setup(p1) and
               IMP.core.RigidMember(p0).get_rigid_body() == IMP.core.RigidMember(p1).get_rigid_body()): continue

            # get distance in the template
            dist0=IMP.core.get_distance(IMP.core.XYZ(p0),IMP.core.XYZ(p1))
            if(dist0<cutoff):
                # define kappa and lower and upper bounds
                kappa=IMP.isd.Scale.setup_particle(IMP.Particle(m),0.0)
                kappa.set_lower(0.0)
                kappa.set_upper(0.0)
                kappa.set_is_optimized(kappa.get_nuisance_key(),False)
                lnar=IMP.isd.LognormalAmbiguousRestraint(p0,p1,kappa,sigmaG)
                lnar.add_contribution(dist0,omega)
                rslist.append(lnar)

    return omega,rslist

def get_kink_restraint(m, rbpairs, kappa):
    rset=IMP.RestraintSet(m, 'Kink_Restraint')
    for rbs in rbpairs:
        h=IMP.core.HarmonicWell([10.0/180.0*math.pi,40.0/180.0*math.pi],kappa)
        kps=IMP.core.RigidBodyAnglePairScore(h)
        pr=IMP.core.PairRestraint(m, kps,IMP.ParticlePair(rbs[0],rbs[1]))
        rset.add_restraint(pr)
    return rset

def get_prior(m, sigmas):
    rslist=[]
    for sigma in sigmas:
        rslist.append(IMP.isd.JeffreysRestraint(m, sigmas[sigma]))
    return rslist

def get_rb_movers(m,rblist,rmax,tmax):
    mvs=[]
    for rb in rblist:
        rblist[rb].set_coordinates_are_optimized(True)
        mv= IMP.core.RigidBodyMover(m, rblist[rb], rmax, tmax)
        mvs.append(mv)
    return mvs

def get_ball_movers(m,ps,tmax):
    mvs=[]
    for p in ps:
        mv= IMP.core.BallMover(m, [p], tmax)
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

def get_crosslink_restraint(m, protein_copies, filename):
    rset=IMP.RestraintSet(m, 'Cysteine_Crosslink')
    # dictionaries
    expdict={}
    epsilons={}

    # beta
    beta=IMP.isd.Scale.setup_particle(IMP.Particle(m),BETA_)
    beta.set_lower(BETA_)
    beta.set_upper(BETA_)
    beta.set_is_optimized(beta.get_nuisance_key(),False)

    # sigma
    sigma=IMP.isd.Scale.setup_particle(IMP.Particle(m),1.0)
    sigma.set_lower(0.0)
    sigma.set_upper(1.0)
    sigma.set_is_optimized(sigma.get_nuisance_key(),False)

    # upper boundary dictionary
    upperb={}
    upperb["16-41"]=1.0-0.628809
    upperb["42-65"]=1.0-0.922462
    upperb["185-226"]=1.0-0.693342

    # epsilon particles
    for id in ["16-41","42-65","185-226"]:
        epsilonscale=IMP.isd.Scale.setup_particle(IMP.Particle(m),0.01)
        epsilonscale.set_lower(0.01)
        epsilonscale.set_upper(upperb[id])
        epsilonscale.set_is_optimized(epsilonscale.get_nuisance_key(),True)
        epsilons[id]=epsilonscale

    # population particle
    pw=IMP.Particle(m)
    weight=IMP.isd.Weight.setup_particle(pw)
    weight.set_weights_are_optimized(True)

    # create grids needed by CrossLinkData
    dist_grid=get_grid(0.0,25.0, NGRID_, False)
    omega1_grid=get_log_grid(1.0, 1000.0, 50)
    sigma_grid=[1.0]
    # read PMF from file
    xlpot=open("../data/cysteine_FES.dat")
    pot_x_grid=[]
    pot_value_grid=[]
    for line in xlpot:
        t=line.split()
        if t[0][0]!="#":
            x = float(t[0])*10.0
            pot_x_grid.append(x)
            pot_value_grid.append(float(t[1])/4.0/math.pi/x/x)
    # CrossLinkMSData
    crossdata=IMP.isd.CrossLinkData(dist_grid,omega1_grid,sigma_grid,pot_x_grid,pot_value_grid,10.0,20.0)

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
        ccldata=IMP.isd.CysteineCrossLinkData(fexp,fmod_grid,omega2_grid,[BETA_])

        # get upperb id
        id=""
        if( 16 <= resid1 <= 41):   id="16-41"
        if( 42 <= resid1 <= 65):   id="42-65"
        if(185 <= resid1 <=226): id="185-226"

        ccl=IMP.isd.CysteineCrossLinkRestraint(m, beta, sigma, epsilons[id], pw,
                                               crossdata, ccldata)

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

def get_Ez_potential(m, protein,boundaries):
    rset=IMP.RestraintSet(m, 'Ez_potential')
    ps=[]
    for boundary in boundaries:
        s=IMP.atom.Selection(protein, residue_indexes=range(*boundary))
        ps.extend(s.get_selected_particle_indexes())
    ez=IMP.atom.EzRestraint(m, ps)
    rset.add_restraint(ez)
    return rset

def initialize_coordinates(m, protein):
    rb=IMP.atom.create_rigid_body(protein)
    rbcoord=rb.get_coordinates()
    rot=IMP.algebra.get_identity_rotation_3d()
    tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
    trans=tmptrans.get_inverse()
    IMP.core.transform(rb,trans)
    rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(IMP.algebra.Transformation3D
       (IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(1,0,0),math.pi/2.0) * IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0,1,0),-math.pi/2.0),
        IMP.algebra.Vector3D(0,0,+14.0))))
    IMP.core.RigidBody.teardown_particle(rb)
    m.remove_particle(rb)

def get_layer_restraint(m, protein,resid,zrange,kappa):
    rset=IMP.RestraintSet(m, 'Layer_restraint')
    lsc=IMP.container.ListSingletonContainer(m)
    s=IMP.atom.Selection(protein, residue_index=resid, atom_type=IMP.atom.AT_CA)
    lsc.add(s.get_selected_particle_indexes())
    hw=IMP.core.HarmonicWell(zrange,kappa)
    asc=IMP.core.AttributeSingletonScore(hw,IMP.FloatKey("z"))
    sr=IMP.container.SingletonsRestraint(asc, lsc)
    rset.add_restraint(sr)
    return rset
