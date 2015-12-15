#!/usr/bin/env python

import sys
import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core
import IMP.rmf
import RMF
import math
sys.path.append('..')
from phoq_util import *

# PARAMETERS
W_BIAS_STRIDE=10000
# max step for Monte Carlo moves
MAXOMEGA_=0.2
MAXTRANS_=0.15
MAXROT_=0.03
MAXEPSILON_=0.01
MAXWEIGHT_=0.1

NITER_=700000
# Run fewer iterations when testing
if '--test' in sys.argv:
    NITER_=500
NOPT_=500
HOMOLOGY_CUTOFF_=8.0
ELASTIC_CUTOFF_=6.0
W_STRIDE=50
KAPPA_=100.0
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

# MODEL
m = IMP.Model()

# set up replica exchange
rem=ReplicaExchange()
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

# read initial model
PHOQ=[]
#filelist=["model_A.pdb","model_B.pdb"]
filelist=["../data/PHOQ_PERI_initial_model.pdb"] * 3
for i in range(0,NCOPIES_):
    PHOQ.append(load_pdb(m, "PHOQ copy "+str(i),filelist[i]))

# initialize system position
for i in range(0,NCOPIES_): initialize_coordinates(m, PHOQ[i])

# create rigid bodies
rblist={}
fixres=[(13,41),(45,184),(194,205),(208,217),(220,233),(245,265)]
for i in range(0,NCOPIES_):
    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
        for res in fixres:
            rb=create_rigid_body(m, chain,res[0],res[1])
            rblist[IMP.atom.Chain(chain).get_id()+"_"+str(res[0])+"-"+str(res[1])+"::"+str(i)]=rb

# Create GLOBAL restraints set
global_rset={}

# Excluded Volume
for i in range(0,NCOPIES_):
    global_rset['Excluded_Volume::'+str(i)]=get_excluded_volume(m, PHOQ[i],1.0)

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
    global_rset['CAForceField::'+str(i)]=IMP.RestraintSet(m)
    excvol_filter['Excluded_Volume::'+str(i)]=IMP.container.ListPairContainer(m)
    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
        for res in floppyres:
            (rstlist,pairslist)=get_CA_force_field(m, chain,res,dihe_dict,ang_dict,do_mix)
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
sigmaG=IMP.isd.Scale.setup_particle(IMP.Particle(m),SIGMAG_)
sigmaG.set_lower(SIGMAG_)
sigmaG.set_upper(SIGMAG_)
sigmaG.set_is_optimized(sigmaG.get_nuisance_key(),False)
ISD_particles.append(sigmaG) # add to list of ISD particles for rmf I/O

# Elastic Network
#elastic_omegas_copies=[]
#elasticres=[(37,49),(181,198),(213,225)]
#elasticres=[(37,49),(181,198),(203,213),(213,225)]
#elasticres=[(37,49),(203,213),(213,225)]

#for i in range(0,NCOPIES_):
#    elastic_omegas={}
#    global_rset['Elastic_Network::'+str(i)]=IMP.RestraintSet(m)
#    for chain in IMP.atom.get_by_type(PHOQ[i], IMP.atom.CHAIN_TYPE):
#        for res in elasticres:
#            (omega,rs)=get_elastic_restraint(m, chain, [res], sigmaG, ELASTIC_CUTOFF_)
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
#template1=load_pdb(m, "template1","../data/3ZRX.pdb")
# 2) 2Y20_A_B
template2=load_pdb(m, "template2","../data/2Y20_A_B.pdb")

# PERIPLASMIC
#peri_align_copies=[]
#peri_kappas_copies=[]
#peri_omegas_copies=[]
#peri_dist_dict_copies=[]
#template3=load_pdb(m, "template3","../data/3BQ8.pdb")

for i in range(0,NCOPIES_):
    hamp_align={}
    #hamp_align['3ZRX']=read_alignments("../data/PHOQ-3ZRX.align",PHOQ[i],template1)
    hamp_align['2Y20']=read_alignments("../data/PHOQ-2Y20.align",PHOQ[i],template2)
    hamp_align_copies.append(hamp_align)

    (hamp_kappas,hamp_omegas,rset,hamp_dist_dict)=get_homology_restraint(m, PHOQ[i],hamp_align,sigmaG,HOMOLOGY_CUTOFF_)
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

    #(peri_kappas,peri_omegas,rset,peri_dist_dict)=get_homology_restraint(m, PHOQ[i],peri_align,sigmaG,HOMOLOGY_CUTOFF_)
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
(global_rset['Cross-link'],crosslink_expdict,crosslink_beta,crosslink_epsilons,crosslink_pw)=get_crosslink_restraint(m, PHOQ,"expcrosslink_90_percent.dat")
# add to list of ISD particles for rmf I/O
ISD_particles.append(crosslink_beta)
ISD_particles.append(crosslink_pw)
for key in crosslink_epsilons:
    ISD_particles.append(crosslink_epsilons[key])

# Jeffrey Prior for omegas
for i in range(0,NCOPIES_):
    global_rset['Prior::'+str(i)]=IMP.RestraintSet(m)
    global_rset['Prior::'+str(i)].add_restraints(get_prior(m, hamp_omegas_copies[i]))
    #global_rset['Prior::'+str(i)].add_restraints(get_prior(m, peri_omegas_copies[i]))
    #global_rset['Prior::'+str(i)].add_restraints(get_prior(m, elastic_omegas_copies[i]))

# Boundaries for Weight
global_rset['Weight Boundaries']=IMP.isd.WeightRestraint(crosslink_pw,0.2,0.8,100000.)

# Ez potential
#TM_regions=[[(13,45)],[(194,221)]]
TM_regions=[(17,37),(195,215)]
for i in range(0,NCOPIES_):
    global_rset['Ez_potential::'+str(i)]=get_Ez_potential(m, PHOQ[i],TM_regions)

# Layer restraint
N_terminus=17
layer=(-17.0,-13.0)
for i in range(0,NCOPIES_):
    global_rset['Layer_restraint::'+str(i)]=get_layer_restraint(m, PHOQ[i],N_terminus,layer,100.0*KAPPA_)

# kink restraint
#for i in range(0,NCOPIES_):
#    hpairs=[]
#    hpairs.append([rblist["A_194-205::"+str(i)],rblist["A_208-217::"+str(i)]])
#    hpairs.append([rblist["B_194-205::"+str(i)],rblist["B_208-217::"+str(i)]])
#    global_rset['Kink::'+str(i)]=get_kink_restraint(m, hpairs, 100.0*KAPPA_)

# make scoring function
sf = IMP.core.RestraintsScoringFunction(list(global_rset.values()))

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
        if(IMP.core.RigidMember.get_is_setup(p)==False):
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
mvs.append(IMP.isd.WeightMover(crosslink_pw,MAXWEIGHT_))

# SerialMover
smv=IMP.core.SerialMover(mvs)

# Monte Carlo with Wte
#mc=IMP.membrane.MonteCarloWithWte(m,EMIN_,EMAX_,SIGMA_,GAMMA_,W0_)
mc=IMP.core.MonteCarlo(m)
mc.set_scoring_function(sf)
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
my_kc=rh.get_category("my data")
score_key=rh.get_key(my_kc, "my score", RMF.float_tag)
cross_key=rh.get_key(my_kc, "my cross-link score", RMF.float_tag)
#bias_key=rh.add_float_key(my_kc,"my bias",True)
index_key=rh.get_key(my_kc, "my index", RMF.int_tag)

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
        myscore=sf.evaluate(False)
        #mybias=mc.get_bias(myscore)
        mycross=global_rset['Cross-link'].evaluate(False)

        # get weights
        ww=IMP.isd.Weight(crosslink_pw).get_weights()

        # prepare printout
        s0="%12d " % (istep)
        s00="%12d " % (myindex)
        #s000="%12.6f " % (mybias)
        s0000="%12.6f " % (mc.get_kt())
        s1="%12.6f " % (mc.get_number_of_accepted_steps()/float(NOPT_))
        s2=' '.join(["%5s %12.6f " % (kkey,global_rset[kkey].evaluate(False)) for kkey in global_rset])
        #s6=' '.join(["%5s %12.6f " % (kkey,sigmas[kkey].get_scale()) for kkey in sigmas])
        #s7=' '.join(["%5s %12.6f " % (kkey,betas[kkey].get_scale()) for kkey in betas])
        s9=' '.join(["%5s %12.6f " % (IMP.isd.CysteineCrossLinkRestraint.get_from(rst).get_name(),IMP.isd.CysteineCrossLinkRestraint.get_from(rst).get_model_frequency()) for rst in global_rset['Cross-link'].get_restraints()])
        s10=' '.join(["%5s %12.6f " % (kkey,crosslink_expdict[kkey]) for kkey in crosslink_expdict])
        s11=' '.join(["%5s %12.6f " % (IMP.isd.CysteineCrossLinkRestraint.get_from(rst).get_name(),IMP.isd.CysteineCrossLinkRestraint.get_from(rst).get_standard_error()) for rst in global_rset['Cross-link'].get_restraints()])
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

        # and other useful data
        rh.get_root_node().set_value(score_key,myscore)
        rh.get_root_node().set_value(cross_key,mycross)
        #rh.get_root_node().set_value(bias_key,mybias)
        rh.get_root_node().set_value(index_key,myindex)

        # print all information of the current frame to rmf
        IMP.rmf.save_frame(rh)

    # time for an exchange
    score=sf.evaluate(False)
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
