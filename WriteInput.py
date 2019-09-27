import yaml
import io
import sys
from math import sqrt, cosh, sinh
import collections


def load_data(): 
    with open("settings.yaml", 'r') as stream:
        data = yaml.safe_load(stream)
    return(data)

def naik(mass):
    mc=float(mass)
    mtree= mc * ( 1 - 3.0*mc**4 / 80.0 + 23*mc**6/2240 + 1783*mc**8/537600 - 76943*mc**10/23654400 )
    epsilon = ( 4 - sqrt( 4 + 12*mtree/( cosh(mtree)*sinh(mtree) ) )) / (sinh(mtree))**2 - 1
    if epsilon > -0.0001:
        return('0.')
    else:
        return('{0:.4f}'.format(epsilon))

def times(data,cfg):
    t0s = []
    nsrc = int(data['lattice info']['nsrc'])
    dsrc = int(data['lattice info']['nt'])/nsrc
    src_start = (19*(cfg/int(data['lattice info']['nsrcdivider']))) % dsrc
    for i in range(nsrc):
        t0s.append(int(src_start + i*dsrc))
    return(t0s)

def make_preamble_2pt(data,input_file,cfg):
    input_file.write('prompt 0\n')
    input_file.write('nx {0}\n'.format(data['lattice info']['nx']))
    input_file.write('ny {0}\n'.format(data['lattice info']['ny']))
    input_file.write('nz {0}\n'.format(data['lattice info']['nz']))
    input_file.write('nt {0}\n'.format(data['lattice info']['nt']))
    input_file.write('iseed {0}\n'.format(cfg))
    input_file.write('job_id hisq-2pt-{0}\n'.format(cfg))
    return()

def linebreak(string,input_file,num):
    input_file.write('\n')
    input_file.write('#')
    input_file.write('='*num)
    input_file.write('\n')
    input_file.write('#  {0}\n'.format(string))
    input_file.write('#')
    input_file.write('='*num)
    input_file.write('\n')
    input_file.write('\n')
    return()

def make_gauge_field(data,cfg,input_file,i):
    input_file.write('#Description of gauge field\n')
    if i == 0:
        input_file.write('reload_serial {0}/{1}.{2}\n'.format(data['lattice info']['cfg_dir'],data['lattice info']['ens'],cfg))
    else:
        input_file.write('continue\n')
    input_file.write('u0 {0}\n'.format(data['lattice info']['u0']))
    input_file.write('no_gauge_fix\n')
    input_file.write('forget\n')
    input_file.write('staple_weight 0\n')
    input_file.write('ape_iter 0\n')
    input_file.write('coordinate_origin 0 0 0 0\n\n')
    input_file.write('#Chiral condensate and related measurements\n\n')
    input_file.write('number_of_pbp_masses 0\n')
    return()

def make_base_source(source_type,source_number,input_file,data,t0,cfg):
    if source_type == 'vec_prop':
        input_file.write('\n')
        input_file.write('\n')
        input_file.write('#Source {0}\n'.format(source_number))
        input_file.write('vector_propagator_file\n')
        input_file.write('subset full\n')
        input_file.write('ncolor {0}\n'.format(data['lattice info']['ncolor']))
        input_file.write('source_label vp\n')
        input_file.write('forget_source')
    if source_type == 'rcw':
        input_file.write('\n')
        input_file.write('\n')
        input_file.write('#Source {0}\n'.format(source_number))
        input_file.write('random_color_wall\n')
        input_file.write('subset full\n')
        input_file.write('t0 {0}\n'.format(t0))
        input_file.write('ncolor {0}\n'.format(data['lattice info']['ncolor']))
        input_file.write('momentum 0 0 0\n')
        input_file.write('source_label rc\n')
        input_file.write('save_serial_scidac_ks_source ./temp/{0}.{1}_t{2}_wallsrc\n'.format(data['lattice info']['ens'],cfg,t0))
    if source_type == 'vec_field':
        input_file.write('\n')
        input_file.write('#Source {0}\n'.format(source_number))
        input_file.write('vector_field\n')
        input_file.write('subset full\n')
        input_file.write('origin 0 0 0 {0}\n'.format(t0))
        input_file.write('load_source ./temp/{0}.{1}_t{2}_wallsrc\n'.format(data['lattice info']['ens'],cfg,t0))
        input_file.write('ncolor {0}\n'.format(data['lattice info']['ncolor']))
        input_file.write('momentum 0 0 0\n')
        input_file.write('source_label vf\n')
        input_file.write('forget_source\n')
    return()

def make_modified_source(input_file,spin_taste,base_source_no,source_number):
    input_file.write('\n')
    input_file.write('\n')
    input_file.write('#Source {0}\n'.format(source_number))
    input_file.write('source {0}\n'.format(base_source_no))
    input_file.write('spin_taste\n')
    input_file.write('spin_taste {0}\n'.format(spin_taste))
    input_file.write('op_label mod\n')
    input_file.write('forget_source\n')
    return()

def make_parent_set(Props,input_file,data,set_no,source,spin_taste,num_props):
    for i in range(num_props):
        Props['st'].append(spin_taste)
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for parent set {0} {1}\n'.format(set_no,spin_taste))
    input_file.write('max_cg_iterations {0}\n'.format(data['lattice info']['max_cg_iterations']))
    input_file.write('max_cg_restarts {0}\n'.format(data['lattice info']['max_cg_restarts']))
    input_file.write('check yes\n')
    input_file.write('momentum_twist {0} {0} {0}\n'.format(data['parent prop']['twist']))
    input_file.write('time_bc periodic\n')
    input_file.write('precison {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators {0}\n'.format(num_props))
    return()

def make_parent_prop(Props,input_file,data,mass,prop_no):
    Props['mass'].append(mass)
    Props['twist'].append(data['parent prop']['twist'])
    Props['type'].append('parent')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))
    input_file.write('mass {0}\n'.format(mass))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(mass)))
    input_file.write('error_for_propagator {0}\n'.format(data['parent prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['parent prop']['rel_error']))
    input_file.write('fresh_ksprop\n')
    input_file.write('forget_ksprop\n')
    input_file.write('\n')    
    return()


def make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_no,prop_no,source,spin_taste,cfg,t0):
    Props['mass'].append(data['daughter prop']['mass'])
    Props['st'].append(spin_taste)
    Props['twist'].append(twist)
    Props['type'].append('daughter')
    if load == True:
        load_directory = data['spectator prop']['load_directory'][data['daughter load']['twists'].index(twist)].format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'])
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for daughter set {0} {1}\n'.format(set_no,spin_taste))
    if load == False:
        input_file.write('max_cg_iterations {0}\n'.format(data['lattice info']['max_cg_iterations']))
        input_file.write('max_cg_restarts {0}\n'.format(data['lattice info']['max_cg_restarts']))
        input_file.write('check yes\n')
    else:
        input_file.write('max_cg_iterations 2\n')
        input_file.write('max_cg_restarts 2\n')
        if check == False:
            input_file.write('check no\n')
        if check == True:
            input_file.write('check yes\n')
    input_file.write('momentum_twist {0} {0} {0}\n'.format(twist))
    input_file.write('time_bc periodic\n')
    input_file.write('precison {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators 1\n')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))                 
    input_file.write('mass {0}\n'.format(data['daughter prop']['mass']))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(data['daughter prop']['mass'])))
    input_file.write('error_for_propagator {0}\n'.format(data['daughter prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['daughter prop']['rel_error']))
    if load == False:
        input_file.write('fresh_ksprop\n')
    else:
        input_file.write('reload_serial_ksprop {0}\n'.format(load_directory))
    if save == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist,spin_taste) )
    else:
        input_file.write('save_serial_scidac_ksprop props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['daughter prop']['mass'],twist,t0))
    input_file.write('\n')    
    return()

def make_spectator_set_prop(Props,input_file,data,set_no,prop_no,source,cfg,t0):
    Props['st'].append('G5-G5')
    Props['mass'].append(data['spectator prop']['mass'])
    Props['twist'].append(data['spectator prop']['twist'])
    Props['type'].append('daughter')  # daughter because don't use in meson unless duplicated in daughter
    load = data['spectator prop']['load']
    save = data['spectator prop']['save']
    check = data['spectator prop']['check']
    if load == True:
        load_directory = data['spectator prop']['load_directory'].format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'])
    spin_taste = 'G5-G5'
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for spectator set {0} {1}\n'.format(set_no,spin_taste))
    if load == False:
        input_file.write('max_cg_iterations {0}\n'.format(data['lattice info']['max_cg_iterations']))
        input_file.write('max_cg_restarts {0}\n'.format(data['lattice info']['max_cg_restarts']))
        input_file.write('check yes\n')
    else:
        input_file.write('max_cg_iterations 2\n')
        input_file.write('max_cg_restarts 2\n')
        if check == False:
            input_file.write('check no\n')
        if check == True:
            input_file.write('check yes\n')   
    input_file.write('momentum_twist {0} {0} {0}\n'.format(data['spectator prop']['twist']))
    input_file.write('time_bc periodic\n')
    input_file.write('precison {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators 1\n')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))                 
    input_file.write('mass {0}\n'.format(data['spectator prop']['mass']))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(data['spectator prop']['mass'])))
    input_file.write('error_for_propagator {0}\n'.format(data['spectator prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['spectator prop']['rel_error']))
    if load == False:
        input_file.write('fresh_ksprop\n')
    else:
        input_file.write('reload_serial_ksprop {0}\n'.format(load_directory))
    if save == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'],'G5-G5') )
    else:
        input_file.write('save_serial_scidac_ksprop props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['spectator prop']['mass'],data['spectator prop']['twist'],t0))
    input_file.write('\n')    
    return()
   
def make_mesons(data,filename,input_file,mass1,mass2,prop1,prop2,twist,t0,spin_taste):   
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Meson masses {0} {1} twist {2} {3}\n'.format(mass1,mass2,twist,spin_taste))
    input_file.write('pair {0} {1}\n'.format(prop1,prop2))
    input_file.write('spectrun_request meson\n')
    input_file.write('save_corr_fnal {0}\n'.format(filename))
    input_file.write('r_offset 0 0 0 {0}\n'.format(t0))
    input_file.write('number_of_correlators 1\n')
    input_file.write('correlator ps p000 1 / {0} {1} 0 0 0 E E E\n'.format(data['lattice info']['wpnorm'],spin_taste))
    return()
    

def make_quarks(num_quarks,input_file):
    input_file.write('\n')
    for i in range(num_quarks):
        input_file.write('\n')
        input_file.write('propagator {0}\n'.format(i))
        input_file.write('identity\n')
        input_file.write('op_label id\n')
        input_file.write('forget_ksprop\n')
    return()

def remove_duplicates(array):
    new_array = []
    new_array.append(array[0])
    for element in array:
        if element not in new_array:
            new_array.append(element)
    return(new_array)


def sources(data):
    source0 = 'rcw'
    if 'daughter load' in data and data['daughter load']['check'] == False:
        source0 = 'vec_prop'
    if data['spectator prop']['load'] == True and data['spectator prop']['check'] == False:
        source0 = 'vec_prop'
    modified_source = []
    for element in data['parent prop']['spin_taste']:
        if element != 'G5-G5':
            modified_source.append(element)
    for element in data['daughter prop']['spin_taste']:
        if element != 'G5-G5':
            modified_source.append(element)
    for element in data['spectator prop']['spin_taste']:
        if element != 'G5-G5':
            modified_source.append(element)
    modified_sources = remove_duplicates(modified_source)
    if source0 == 'rcw':
       no_base_sources = 2 
    elif source0 == 'vec_prop':
       no_base_sources = 3
    return(no_base_sources,source0,modified_sources)


def no_sets_mesons(data):
    if 'daughter existing' in data:
        take_off_existing = len(data['daughter existing']['twists'])
        if data['spectator prop']['twist'] in data['daughter existing']['twists'] and data['daughter existing']['spin_taste'][data['daughter existing']['twists'].index(data['spectator prop']['twist'])] == 'G5-G5':
            take_off_existing = take_off_existing - 1
    else:
        take_off_existing = 0
    if data['parent prop']['multimass'] == True:
        if data['spectator prop']['same'] == True:
            set_no = len(data['parent prop']['spin_taste']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing 
        elif data['spectator prop']['same'] == False:
            set_no = len(data['parent prop']['spin_taste']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing + 1
    if data['parent prop']['multimass'] == False:
        if data['spectator prop']['same'] == True:
            set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing 
        elif data['spectator prop']['same'] == False:
            set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing + 1
    meson_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
    if data['spectator prop']['same'] == True:
        prop_no = meson_no
    elif data['spectator prop']['same'] == False:
        prop_no = meson_no + 1 
    return(set_no,meson_no,prop_no)


def main_2pts(argv):
    cfg = int(argv[0])
    input_file = open('./input-2pt/milc_2pt_{0}.in'.format(cfg), 'w+')
    data = load_data()
    t0s = times(data,cfg)
    no_base_sources,source0, modified_sources = sources(data)
    make_preamble_2pt(data,input_file,cfg)
    for i,t0 in enumerate(t0s):
        set_no,meson_no,prop_no = no_sets_mesons(data)
        Props = collections.OrderedDict()
        Props['st'] = []
        Props['mass'] = []
        Props['twist'] = []
        Props['type'] = []
        linebreak('Source time {0}'.format(t0) ,input_file,80)
        make_gauge_field(data,cfg,input_file,i)
        linebreak('Description of base sources',input_file,40)
        ################### BASE SOURCES #####################################
        input_file.write('number_of_base_sources {0}\n'.format(no_base_sources))
        if source0 == 'vec_prop':
            make_base_source('vec_prop',0,input_file,data,t0,cfg)
            make_base_source('rcw',1,input_file,data,t0,cfg)
            make_base_source('vec_field',2,input_file,data,t0,cfg)
        elif source0 == 'rcw': 
            make_base_source('rcw',0,input_file,data,t0,cfg)
            make_base_source('vec_field',1,input_file,data,t0,cfg)
        
        linebreak('Description of modified sources',input_file,40)
        ################### MODIFIED SOURCES #################################
        input_file.write('number_of_modified_sources {0}\n'.format(len(modified_sources)))
        if source0 == 'vec_prop':
            n = 3
            for i in range(len(modified_sources)):
                make_modified_source(input_file,modified_sources[i],2,n)
                n += 1
        elif source0 == 'rcw':
            n = 2
            for i in range(len(modified_sources)):
                make_modified_source(input_file,modified_sources[i],1,n)
                n += 1
                
        linebreak('Description of propagators',input_file,40)
        ################### PROPAGATORS ######################################
        pr_num = 0
        set_num = 0
        input_file.write('number_of_sets {0}\n'.format(set_no))
        #~~~~~~~~~~~~~~~~~~ SPECTATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
        if data['spectator prop']['check'] == True and source0 == 'vec_prop':
            make_spectator_set_prop(Props,input_file,data,set_num,pr_num,1,cfg,t0)
            pr_num+=1
            set_num+=1
        else:
            make_spectator_set_prop(Props,input_file,data,set_num,pr_num,0,cfg,t0)
            pr_num+=1
            set_num+=1
        #~~~~~~~~~~~~~~~~~~ DAUGHTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for st in remove_duplicates(data['daughter prop']['spin_taste']):
            for twist in data['daughter prop']['twists']:
                load = False
                if st == 'G5-G5':
                    save = data['daughter prop']['save']
                else:
                    save = False
                check = True
                
                if 'daughter existing' in data and twist in data['daughter existing']['twists'] and st==data['daughter existing']['spin_taste'][data['daughter existing']['twists'].index(twist)]:
                    pass
                elif twist == data['spectator prop']['twist'] and st == 'G5-G5' and data['spectator prop']['mass']==data['daughter prop']['mass']:
                    pass
                elif 'daughter load' in data and twist in data['daughter load']['twists'] and st==data['daughter load']['spin_taste'][data['daughter load']['twists'].index(twist)]:
                    load = True
                    save = False
                    check = data['daughter load']['check']
                    if data['daughter load']['check'] == True and source0 == 'vec_prop':
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,1,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                    else:
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,0,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                elif st =='G5-G5':
                    if source0 == 'vec_prop':
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,1,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                    else:
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,0,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                else:
                    if source0 == 'vec_prop':
                        sourcenum = 3 + modified_sources.index(st)
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sourcenum,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                    else:
                        sourcenum = 2 + modified_sources.index(st)
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sourcenum,st,cfg,t0)
                        pr_num+=1
                        set_num+=1
        #~~~~~~~~~~~~~~~~~~~~~ PARENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if data['parent prop']['multimass']== False:
            for st in remove_duplicates(data['parent prop']['spin_taste']):
                for mass in data['parent prop']['masses']:
                    if st =='G5-G5':
                        if source0 == 'vec_prop':
                            make_parent_set(Props,input_file,data,set_num,1,st,1)
                        else:
                            make_parent_set(Props,input_file,data,set_num,0,st,1)
                    else:
                        if source0 == 'vec_prop':
                            sourcenum = 3 + modified_sources.index(st)
                            make_parent_set(Props,input_file,data,set_num,sourcenum,st,1)
                        else:
                            sourcenum = 2 + modified_sources.index(st)
                            make_parent_set(Props,input_file,data,set_num,sourcenum,st,1)                            
                    make_parent_prop(Props,input_file,data,mass,pr_num)
                    pr_num+=1
                    set_num+=1
        else:
            for st in remove_duplicates(data['parent prop']['spin_taste']):
                if st =='G5-G5':
                    if source0 == 'vec_prop':
                        make_parent_set(Props,input_file,data,set_num,1,st,len(data['parent prop']['masses']))
                    else:
                        make_parent_set(Props,input_file,data,set_num,0,st,len(data['parent prop']['masses']))
                else:
                    if source0 == 'vec_prop':
                        sourcenum = 3 + modified_sources.index(st)
                        make_parent_set(Props,input_file,data,set_num,sourcenum,st,len(data['parent prop']['masses']))
                    else:
                        sourcenum = 2 + modified_sources.index(st)
                        make_parent_set(Props,input_file,data,set_num,sourcenum,st,len(data['parent prop']['masses'])) 
                for mass in data['parent prop']['masses']:
                    make_parent_prop(Props,input_file,data,mass,pr_num)
                    pr_num+=1
                set_num+=1
            
        linebreak('Description of quarks',input_file,40)
        ################### QUARKS ###########################################
        input_file.write('number_of_quarks {0}\n'.format(prop_no))
        make_quarks(prop_no,input_file)
        linebreak('Description of mesons',input_file,40)
        ################### MESONS ###########################################
        if 'daughter existing' in data and data['spectator prop']['twist'] in data['daughter existing']['twists'] and data['daughter existing']['spin_taste'][data['daughter existing']['twists'].index(data['spectator prop']['twist'])] == 'G5-G5':
            start = 1
            meson_no = meson_no -1
        else:
            start = 0
        input_file.write('number_of_mesons {0}\n'.format(meson_no))
        if data['spectator prop']['same'] == True:
            for i in range(start,meson_no+start):
                filename = data['{0} prop'.format(Props['type'][i])]['corr_files'][data['{0} prop'.format(Props['type'][i])]['spin_taste'].index(Props['st'][i])].format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,Props['mass'][0],Props['mass'][i])
                make_mesons(data,filename,input_file,Props['mass'][0],Props['mass'][i],0,i,Props['twist'][i],t0,Props['st'][i])
        if data['spectator prop']['same'] == False:
            for i in range(1,meson_no+1):
                filename = data['{0} prop'.format(Props['type'][i])]['corr_files'][data['{0} prop'.format(Props['type'][i])]['spin_taste'].index(Props['st'][i])].format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,Props['mass'][0],Props['mass'][i])
                make_mesons(data,filename,input_file,Props['mass'][0],Props['mass'][i],0,i,Props['twist'][i],t0,Props['st'][i])
        
        linebreak('Description of baryons',input_file,40)
        ################### BARYONS ##########################################
        input_file.write('number_of_baryons 0\n')
    input_file.close()
    return()


########################################### EXTSRC #####################################################################

def Times(data,cfg,t0):
    Ts = []
    actual_Ts = []
    for i in range(data['three points']['nT']):
        actual_Ts.append(t0 + data['three points']['Tstart'] + i * data['three points']['dT'])
        Ts.append((t0 + data['three points']['Tstart'] + i * data['three points']['dT']) % data['lattice info']['nt'])
    return(Ts,actual_Ts)

def make_preamble_ext(data,input_file,cfg):
    input_file.write('prompt 0\n')
    input_file.write('nx {0}\n'.format(data['lattice info']['nx']))
    input_file.write('ny {0}\n'.format(data['lattice info']['ny']))
    input_file.write('nz {0}\n'.format(data['lattice info']['nz']))
    input_file.write('nt {0}\n'.format(data['lattice info']['nt']))
    input_file.write('\n')
    input_file.write('job_id hisq-extsrc-{0}\n'.format(cfg))
    return()

def make_quark_ext(input_file,data,qnum,cfg,t0,Ts):
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('# Quark {0}\n'.format(qnum))
    input_file.write('quark_type KS\n')
    input_file.write('output_type KS\n')
    if data['spectator prop']['save'] == True:
        input_file.write('reload_serial_ksprop props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['spectator prop']['mass'],data['spectator prop']['twist'],t0))
    else:
        input_file.write('reload_serial_ksprop ./temp/{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'],'G5-G5') )
    input_file.write('ncolor {0}\n'.format(data['lattice info']['ncolor']))
    input_file.write('\n')
    input_file.write('# Smeraing for quark {0}\n'.format(qnum))
    input_file.write('identity\n')
    input_file.write('op_label ext\n')
    input_file.write('sink_gamma {0}\n'.format(data['parent prop']['spin_taste'][qnum]))
    input_file.write('\n')
    input_file.write('r_offset 0 0 0 0\n')
    input_file.write('number_of_time_slices {0}\n'.format(data['three points']['nT']))
    for T in Ts:
        input_file.write('save_serial_scidac_ks_source ./temp/{0}.{1}_t{2}_extsrc_{3}_T{4}_m{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['parent prop']['spin_taste'][qnum],T,data['spectator prop']['mass']))
        input_file.write('t0 {0}\n'.format(T))
    return()

def main_extsrc(argv):
    cfg = int(argv[0])
    input_file = open('./input-extsrc/milc_ext_{0}.in'.format(cfg), 'w+')
    data = load_data()
    t0s = times(data,cfg)
    make_preamble_ext(data,input_file,cfg)
    for t0 in t0s:
        Ts,actual_Ts = Times(data,cfg,t0)
        linebreak('Source time {0}'.format(t0) ,input_file,80)
        input_file.write('number_of_quarks {0}\n'.format(len(data['parent prop']['spin_taste'])))
        for i in range(len(data['parent prop']['spin_taste'])):
            make_quark_ext(input_file,data,i,cfg,t0,Ts)
    return()

############################################# 3 pts ###########################################################


def make_preamble_3pt(data,input_file,cfg):
    input_file.write('prompt 0\n')
    input_file.write('nx {0}\n'.format(data['lattice info']['nx']))
    input_file.write('ny {0}\n'.format(data['lattice info']['ny']))
    input_file.write('nz {0}\n'.format(data['lattice info']['nz']))
    input_file.write('nt {0}\n'.format(data['lattice info']['nt']))
    input_file.write('iseed {0}\n'.format(cfg))
    input_file.write('job_id hisq-3pt-{0}\n'.format(cfg))
    return()



def main_3pts(argv):
    cfg = int(argv[0])
    input_file = open('./input-3pt/milc_3pt_{0}.in'.format(cfg), 'w+')
    data = load_data()
    t0s = times(data,cfg)
    make_preamble_3pt(data,input_file,cfg)
    i = 0
    for t0 in t0s:
        Ts,actual_Ts = Times(data,cfg,t0)
        for j,T in enumerate(Ts):
            linebreak('Source time {0}, meson separation T = {1}'.format(t0,actual_Ts[j]-t0) ,input_file,80)
            make_gauge_field(data,cfg,input_file,i)
            i += 1 
            linebreak('Description of base sources',input_file,40)
    return()









main_2pts(sys.argv[1:])
main_extsrc(sys.argv[1:])
main_3pts(sys.argv[1:])
