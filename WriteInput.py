import yaml
import io
import sys
from math import sqrt, cosh, sinh
import collections
import os

reload_check_iters = 1000

def load_data(): 
    with open("./in/settings.yaml", 'r') as stream:
        data = yaml.safe_load(stream)
    return(data)

data = load_data()

def edit_submit(data): # edits submit script so the input file read and written includes the tag
    cfg = int(sys.argv[1:][0])
    f = open(data['lattice info']['submit'], 'r')
    lines = f.readlines()
    f.close()
    if not os.path.exists('./submit'):
        os.makedirs('./submit')
    g = open('./submit/use_this_submit', 'w+')
    for line in lines:
        if (line.find('milc_in') != -1 or line.find('milc_out') != -1) and line.find('mpirun') ==-1 :
            g.write(line[:line.find('/milc')+1] + data['lattice info']['tag'] + line[line.find('/milc')+1:])
        elif line.find('local ens=') != -1:
            g.write('  local ens={0}'.format(data['lattice info']['ens']))
        elif line.find('#SBATCH --job-name=') != -1:
            if 'parent prop' in data and 'daughter prop' in data and data['lattice info']['justtwopoints'] == False:
                g.write('#SBATCH --job-name={0}{1}{2}-{3}\n'.format(data['lattice info']['tag'],data['parent prop']['name'],data['daughter prop']['name'],cfg))
            else:
                g.write('#SBATCH --job-name={0}2pt-{1}\n'.format(data['lattice info']['tag'],cfg))
        elif line.find('#SBATCH --output=') != -1:
            if 'parent prop' in data and 'daughter prop' in data and data['lattice info']['justtwopoints'] == False:
                g.write('#SBATCH --output=./out/{0}{1}{2}-{3}-%A.out\n'.format(data['lattice info']['tag'],data['parent prop']['name'],data['daughter prop']['name'],cfg))
            else:
                g.write('#SBATCH --output=./out/{0}2pt-{1}-%A.out\n'.format(data['lattice info']['tag'],cfg))
        elif line.find('#SBATCH --error=') != -1:
            if 'parent prop' in data and 'daughter prop' in data and data['lattice info']['justtwopoints'] == False:
                g.write('#SBATCH --error=./out/{0}{1}{2}-{3}-%A.err\n'.format(data['lattice info']['tag'],data['parent prop']['name'],data['daughter prop']['name'],cfg))
            else:
                g.write('#SBATCH --error=./out/{0}2pt-{1}-%A.err\n'.format(data['lattice info']['tag'],cfg))
        
        elif line.find('rm ${temp}/${ens}*${cfg}_*') != -1:
            g.write('  rm {0}temp{1}/{2}{0}ens{1}*{0}cfg{1}_*\n'.format('${','}',data['lattice info']['tag']))
        else:
            g.write(line)
    g.close()
    return()

edit_submit(data)

def make_directories(data):
    if not os.path.exists('./corrs'):
        os.makedirs('./corrs')
    if not os.path.exists('./temp'):
        os.makedirs('./temp')
    if not os.path.exists('./props'):
        os.makedirs('./props')
    if not os.path.exists('./out'):
        os.makedirs('./out')
    if not os.path.exists('./in/input-2pt'):
        os.makedirs('./in/input-2pt')
    if data['lattice info']['justtwopoints'] == False:
        if not os.path.exists('./in/input-extsrc'):
            os.makedirs('./in/input-extsrc')
        if not os.path.exists('./in/input-3pt'):
            os.makedirs('./in/input-3pt')
    if 'parent prop' in data:
        for st in data['parent prop']['spin_taste']:
            if not os.path.exists('./corrs/{0}_{1}_tw{2}'.format(data['parent prop']['name'],st,data['parent prop']['twist'])):
                os.makedirs('./corrs/{0}_{1}_tw{2}'.format(data['parent prop']['name'],st,data['parent prop']['twist']))
                
    if 'daughter prop' in data:
        for st in data['daughter prop']['spin_taste']:
            for twist in data['daughter prop']['twists']:
                if not os.path.exists('./corrs/{0}_{1}_tw{2}'.format(data['daughter prop']['name'],st,twist)):
                    os.makedirs('./corrs/{0}_{1}_tw{2}'.format(data['daughter prop']['name'],st,twist))
    if data['lattice info']['justtwopoints'] == False:
        for k in range(len(data['three points']['p J d'])):
            for twist in data['daughter prop']['twists']:
                if not os.path.exists('./corrs/current-{0}_tw{1}'.format(data['three points']['label'][k],twist)):
                    os.makedirs('./corrs/current-{0}_tw{1}'.format(data['three points']['label'][k],twist))
                    
        
    return()
make_directories(data)

def naik(mass,quark):
    mc=float(mass)
    mtree= mc * ( 1 - 3.0*mc**4 / 80.0 + 23*mc**6/2240 + 1783*mc**8/537600 - 76943*mc**10/23654400 )
    epsilon = ( 4 - sqrt( 4 + 12*mtree/( cosh(mtree)*sinh(mtree) ) )) / (sinh(mtree))**2 - 1
    if quark == 'strange' or quark == 'light':
        return('0.')
    else:
        return('{0:.4f}'.format(epsilon))

def times(data,cfg):
    t0s = []
    nsrc = int(data['lattice info']['nsrc'])
    dsrc = int(data['lattice info']['nt'])/nsrc
    src_start = (19*(cfg/int(data['lattice info']['nsrcdivider']))) % dsrc
    if data['lattice info']['allsources'] ==  True:
        for i in range(nsrc):
            t0s.append(int(src_start + i*dsrc))
    else:
        for i in data['lattice info']['allsources']:
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

def make_base_source(source_type,source_number,input_file,data,t0,cfg,filename):
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
        input_file.write('save_serial_scidac_ks_source {0}\n'.format(filename))
    if source_type == 'vec_field':
        input_file.write('\n')
        input_file.write('\n')
        input_file.write('#Source {0}\n'.format(source_number))
        input_file.write('vector_field\n')
        input_file.write('subset full\n')
        input_file.write('origin 0 0 0 {0}\n'.format(t0))
        input_file.write('load_source {0}\n'.format(filename))
        input_file.write('ncolor {0}\n'.format(data['lattice info']['ncolor']))
        input_file.write('momentum 0 0 0\n')
        input_file.write('source_label vf\n')
        input_file.write('forget_source\n')
    return()

def make_modified_source(input_file,spin_taste,source_number,sources):
    input_file.write('\n')
    input_file.write('\n')
    input_file.write('#Source {0}\n'.format(source_number))
    input_file.write('source {0}\n'.format(sources.index('vec_field')))
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
    input_file.write('precision {0}\n'.format(data['lattice info']['precision']))
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
    input_file.write('naik_term_epsilon {0}\n'.format(naik(mass,data['parent prop']['quark'])))
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
        load_directory = data['daughter load']['load_directory'][data['daughter load']['twists'].index(twist)].format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist)
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for daughter set {0} {1}\n'.format(set_no,spin_taste))
    if load == False:
        input_file.write('max_cg_iterations {0}\n'.format(data['lattice info']['max_cg_iterations']))
        input_file.write('max_cg_restarts {0}\n'.format(data['lattice info']['max_cg_restarts']))
        input_file.write('check yes\n')
    else:
        input_file.write('max_cg_iterations {0}\n'.format(reload_check_iters))
        input_file.write('max_cg_restarts 3\n')
        if check == False:
            input_file.write('check no\n')
        if check == True:
            input_file.write('check yes\n')
    input_file.write('momentum_twist {0} {0} {0}\n'.format(twist))
    input_file.write('time_bc periodic\n')
    input_file.write('precision {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators 1\n')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))                 
    input_file.write('mass {0}\n'.format(data['daughter prop']['mass']))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(data['daughter prop']['mass'],data['daughter prop']['quark'])))
    input_file.write('error_for_propagator {0}\n'.format(data['daughter prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['daughter prop']['rel_error']))
    if load == False:
        input_file.write('fresh_ksprop\n')
    else:
        input_file.write('reload_serial_ksprop {0}\n'.format(load_directory))
    if save == False and data['lattice info']['justtwopoints'] == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist,spin_taste,data['lattice info']['tag']) )
    elif spin_taste != 'G5-G5' and data['lattice info']['justtwopoints'] == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist,spin_taste,data['lattice info']['tag']) )
    elif save == True:
        input_file.write('save_serial_scidac_ksprop ./props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['daughter prop']['mass'],twist,t0))
    else:
        input_file.write('forget_ksprop\n')
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
        load_directory = data['spectator prop']['load_directory'].format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'])
    spin_taste = 'G5-G5'
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for spectator set {0} {1}\n'.format(set_no,spin_taste))
    if load == False:
        input_file.write('max_cg_iterations {0}\n'.format(data['lattice info']['max_cg_iterations']))
        input_file.write('max_cg_restarts {0}\n'.format(data['lattice info']['max_cg_restarts']))
        input_file.write('check yes\n')
    else:
        input_file.write('max_cg_iterations {0}\n'.format(reload_check_iters))
        input_file.write('max_cg_restarts 3\n')
        if check == False:
            input_file.write('check no\n')
        if check == True:
            input_file.write('check yes\n')   
    input_file.write('momentum_twist {0} {0} {0}\n'.format(data['spectator prop']['twist']))
    input_file.write('time_bc periodic\n')
    input_file.write('precision {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators 1\n')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))                 
    input_file.write('mass {0}\n'.format(data['spectator prop']['mass']))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(data['spectator prop']['mass'],data['spectator prop']['quark'])))
    input_file.write('error_for_propagator {0}\n'.format(data['spectator prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['spectator prop']['rel_error']))
    if load == False:
        input_file.write('fresh_ksprop\n')
    else:
        input_file.write('reload_serial_ksprop {0}\n'.format(load_directory))
    if save == False and data['lattice info']['justtwopoints'] == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'],'G5-G5',data['lattice info']['tag']) )
    elif load == True and data['lattice info']['justtwopoints'] == False:
        input_file.write('save_serial_scidac_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'],'G5-G5',data['lattice info']['tag']) )
    elif save == True:
        input_file.write('save_serial_scidac_ksprop ./props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['spectator prop']['mass'],data['spectator prop']['twist'],t0))
        
    else:
       input_file.write('forget_ksprop\n')
    input_file.write('\n')    
    return()
   
def make_mesons(data,filename,input_file,mass1,mass2,prop1,prop2,twist,t0,spin_taste):   
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Meson masses {0} {1} twist {2} {3}\n'.format(mass1,mass2,twist,spin_taste))
    input_file.write('pair {0} {1}\n'.format(prop1,prop2))
    input_file.write('spectrum_request meson\n')
    input_file.write('save_corr_fnal ./corrs/{0}\n'.format(filename))
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
    if len(array) > 0:
        new_array.append(array[0])
        for element in array:
            if element not in new_array:
                new_array.append(element)
    return(new_array)

def num_G5sets_gen(data):
    num_gen = 0
    if 'parent prop' in data and data['parent prop']['multimass'] == False and 'G5-G5' in data['parent prop']['spin_taste']:
        num_gen += len(data['parent prop']['masses'])
    if 'parent prop' in data and data['parent prop']['multimass'] == True and 'G5-G5' in data['parent prop']['spin_taste']:
        num_gen += 1
    if 'daughter prop' in data and 'G5-G5' in data['daughter prop']['spin_taste']:
        num_gen += len(data['daughter prop']['twists'])
    if 'daughter existing' in data:
        for ste in data['daughter existing']['spin_taste']:
            if ste == 'G5-G5':
                num_gen -= 1
    if 'daughter load' in data:
        for stl in data['daughter load']['spin_taste']:
            if stl == 'G5-G5':
                num_gen -= 1
    if 'daughter prop' in data and data['spectator prop']['same'] == False and data['spectator prop']['load'] == False:
        num_gen += 1
    if 'daughter prop' not in data and data['spectator prop']['load'] == False:
        num_gen += 1
    return(num_gen)


def which_sources(data):
    num_gen = num_G5sets_gen(data)
    sources = []
    if 'daughter load' in data and data['daughter load']['check'] == False:
        sources.append('vec_prop')
    elif data['spectator prop']['load'] == True and data['spectator prop']['check'] == False:
        sources.append('vec_prop')    
    modified_source = []
    if 'parent prop' in data:
        for element in data['parent prop']['spin_taste']:
            if element != 'G5-G5':
                modified_source.append(element)
    if 'daughter prop' in data:
        for element in data['daughter prop']['spin_taste']:
            if element != 'G5-G5':
                modified_source.append(element)
    modified_sources = remove_duplicates(modified_source)
    if num_gen == 1 and len(modified_sources) == 0:
        sources.append('rcw')
    elif num_gen == 0 and len(modified_sources) == 0:
        pass
    else:
        sources.append('rcw')
        sources.append('vec_field')
    return(sources,modified_sources)


def no_sets_mesons(data):
    if 'daughter existing' in data:
        take_off_existing = len(data['daughter existing']['twists'])
        if data['spectator prop']['twist'] in data['daughter existing']['twists'] and data['daughter existing']['spin_taste'][data['daughter existing']['twists'].index(data['spectator prop']['twist'])] == 'G5-G5' and data['spectator prop']['same'] == True:
            take_off_existing = take_off_existing - 1
    else:
        take_off_existing = 0
    if 'parent prop' in data:
        if data['parent prop']['multimass'] == True:
            if data['spectator prop']['same'] == True:
                if 'daughter prop' in data: 
                    set_no = len(data['parent prop']['spin_taste']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
                else:
                    set_no = len(data['parent prop']['spin_taste']) + 1
            else:
                if 'daughter prop' in data:
                    set_no = len(data['parent prop']['spin_taste']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing + 1
                else:
                    set_no = len(data['parent prop']['spin_taste']) + 1
        if data['parent prop']['multimass'] == False:
            if data['spectator prop']['same'] == True:
                if 'daughter prop' in data: 
                    set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
                else:
                    set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses'])  
            else:
                if 'daughter prop' in data: 
                    set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing + 1
                else:
                    set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + 1
                    
    else:
        set_no = len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
                

    if 'parent prop' in data and 'daughter prop' in data:
        meson_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
    elif 'daughter prop' in data:
        meson_no =  len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) - take_off_existing
    elif 'parent prop' in data:
        meson_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses'])  - take_off_existing
       
    if data['spectator prop']['same'] == True:
        prop_no = meson_no
    elif data['spectator prop']['same'] == False:
        prop_no = meson_no + 1 
    return(set_no,meson_no,prop_no)


def main_2pts(argv):
    data = load_data()
    cfg = int(argv[0])
    input_file = open('./in/input-2pt/{0}milc_2pt_{1}.in'.format(data['lattice info']['tag'],cfg), 'w+')
    t0s = times(data,cfg)
    sources, modified_sources = which_sources(data)
    make_preamble_2pt(data,input_file,cfg)
    for i,t0 in enumerate(t0s):
        set_no,meson_no,prop_no = no_sets_mesons(data)
        Props = collections.OrderedDict()
        Props['st'] = []
        Props['mass'] = []
        Props['twist'] = []
        Props['type'] = []
        linebreak('Source time {0} ({1} of {2})'.format(t0,i+1,len(t0s)) ,input_file,80)
        make_gauge_field(data,cfg,input_file,i)
        linebreak('Description of base sources',input_file,40)
        ################### BASE SOURCES #####################################            
        input_file.write('number_of_base_sources {0}\n'.format(len(sources)))
        fname = './temp/{3}{0}.{1}_t{2}_wallsrc'.format(data['lattice info']['ens'],cfg,t0,data['lattice info']['tag'])
        for base_source_no,source in enumerate(sources):
            make_base_source(source,base_source_no,input_file,data,t0,cfg,fname)        
        linebreak('Description of modified sources',input_file,40)
        ################### MODIFIED SOURCES #################################
        input_file.write('number_of_modified_sources {0}\n'.format(len(modified_sources)))
        for mod_source_no,source in enumerate(modified_sources):
            make_modified_source(input_file,source,len(sources)+mod_source_no,sources) 
        linebreak('Description of propagators',input_file,40)
        ################### PROPAGATORS ######################################
        pr_num = 0
        set_num = 0
        gen_prop_no = 0
        input_file.write('number_of_sets {0}\n'.format(set_no))
        #~~~~~~~~~~~~~~~~~~ SPECTATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
        if data['spectator prop']['load'] == False:
            make_spectator_set_prop(Props,input_file,data,set_num,pr_num,sources.index('rcw'),cfg,t0)
            pr_num+=1
            set_num+=1
            gen_prop_no = 1
        elif data['spectator prop']['check'] == True:
            make_spectator_set_prop(Props,input_file,data,set_num,pr_num,sources.index('rcw'),cfg,t0)
            pr_num+=1
            set_num+=1
        else:
            make_spectator_set_prop(Props,input_file,data,set_num,pr_num,sources.index('vec_prop'),cfg,t0)
            pr_num+=1
            set_num+=1
        #~~~~~~~~~~~~~~~~~~ DAUGHTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if 'daughter prop' in data:
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
                    elif twist == data['spectator prop']['twist'] and st == 'G5-G5' and data['spectator prop']['mass']==data['daughter prop']['mass'] and data['spectator prop']['same'] == True:
                        pass
                    elif 'daughter load' in data and twist in data['daughter load']['twists'] and st==data['daughter load']['spin_taste'][data['daughter load']['twists'].index(twist)]:
                        load = True
                        save = False
                        check = data['daughter load']['check']
                        if data['daughter load']['check'] == True:
                            make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sources.index('rcw'),st,cfg,t0)
                            pr_num+=1
                            set_num+=1
                        else:
                            make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sources.index('vec_prop'),st,cfg,t0)
                            pr_num+=1
                            set_num+=1
                    elif st =='G5-G5':
                        if gen_prop_no == 0:
                            make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sources.index('rcw'),st,cfg,t0)
                            pr_num+=1
                            set_num+=1
                            gen_prop_no = 1
                        else:
                            make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,sources.index('vec_field'),st,cfg,t0)
                            pr_num+=1
                            set_num+=1
                    else:
                        make_daughter_set_prop(Props,load,save,check,input_file,data,twist,set_num,pr_num,len(sources)+modified_sources.index(st),st,cfg,t0)
                        pr_num+=1
                        set_num+=1
                    
        #~~~~~~~~~~~~~~~~~~~~~ PARENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if 'parent prop' in data:
            if data['parent prop']['multimass']== False:
                for st in remove_duplicates(data['parent prop']['spin_taste']):
                    for mass in data['parent prop']['masses']:
                        if st =='G5-G5':
                            if gen_prop_no == 0:
                                make_parent_set(Props,input_file,data,set_num,sources.index('rcw'),st,1)
                                gen_prop_no = 1
                            else:
                                make_parent_set(Props,input_file,data,set_num,sources.index('vec_field'),st,1)
                        else:
                            make_parent_set(Props,input_file,data,set_num,len(sources) + modified_sources.index(st),st,1)                   
                        make_parent_prop(Props,input_file,data,mass,pr_num)
                        pr_num+=1
                        set_num+=1
            else:
                for st in remove_duplicates(data['parent prop']['spin_taste']):
                    if st =='G5-G5':
                        if gen_prop_no == 0:
                            make_parent_set(Props,input_file,data,set_num,sources.index('rcw'),st,len(data['parent prop']['masses']))
                            gen_prop_no = 1
                        else:
                            make_parent_set(Props,input_file,data,set_num,sources.index('vec_field'),st,len(data['parent prop']['masses']))
                    else:
                        make_parent_set(Props,input_file,data,set_num,len(sources) + modified_sources.index(st),st,len(data['parent prop']['masses']))
                       
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
                if Props['type'][i] == 'parent':
                    filename = '{8}_{7}_tw{0}/{1}meson-{8}-{7}.{2}.{3}_t{4}_m{5}_m{6}'.format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],Props['mass'][i],Props['st'][i],data['parent prop']['name'])
                elif Props['type'][i] == 'daughter':
                    filename = '{8}_{7}_tw{0}/{1}meson-{8}-{7}.{2}.{3}_t{4}_m{5}_m{6}'.format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],Props['mass'][i],Props['st'][i],data['daughter prop']['name'])
                    
                make_mesons(data,filename,input_file,Props['mass'][0],Props['mass'][i],0,i,Props['twist'][i],t0,Props['st'][i])
        if data['spectator prop']['same'] == False:
            for i in range(1,meson_no+1):
                if Props['type'][i] == 'parent':
                    filename = '{8}_{7}_tw{0}/{1}meson-{8}-{7}.{2}.{3}_t{4}_m{5}_m{6}'.format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],Props['mass'][i],Props['st'][i],data['parent prop']['name'])
                elif Props['type'][i] == 'daughter':
                    filename = '{8}_{7}_tw{0}/{1}meson-{8}-{7}.{2}.{3}_t{4}_m{5}_m{6}'.format(Props['twist'][i],data['lattice info']['tag'],data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],Props['mass'][i],Props['st'][i],data['daughter prop']['name'])
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
    if data['spectator prop']['save'] == True and data['spectator prop']['load'] == False:
        input_file.write('reload_serial_ksprop ./props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['spectator prop']['mass'],data['spectator prop']['twist'],t0))
    else:
        input_file.write('reload_serial_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['spectator prop']['mass'],data['spectator prop']['twist'],'G5-G5',data['lattice info']['tag']) )
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
        input_file.write('save_serial_scidac_ks_source ./temp/{6}{0}.{1}_t{2}_extsrc_{3}_T{4}_m{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['parent prop']['spin_taste'][qnum],T,data['spectator prop']['mass'],data['lattice info']['tag']))
        input_file.write('t0 {0}\n'.format(T))
    return()

def main_extsrc(argv):
    data = load_data()
    cfg = int(argv[0])
    input_file = open('./in/input-extsrc/{0}milc_ext_{1}.in'.format(data['lattice info']['tag'],cfg), 'w+')
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

def make_daughter_set_prop_3pt(Props,save,input_file,data,twist,set_no,prop_no,source,spin_taste,cfg,t0):
    Props['mass'].append(data['daughter prop']['mass'])
    Props['st'].append(spin_taste)
    Props['twist'].append(twist)
    Props['type'].append('daughter')
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('#Parameters for daughter set {0} {1}\n'.format(set_no,spin_taste))
    input_file.write('max_cg_iterations 2\n')
    input_file.write('max_cg_restarts 2\n')
    input_file.write('check no\n')
    input_file.write('momentum_twist {0} {0} {0}\n'.format(twist))
    input_file.write('time_bc periodic\n')
    input_file.write('precision {0}\n'.format(data['lattice info']['precision']))
    input_file.write('\n')
    input_file.write('source {0}\n'.format(source))
    input_file.write('number_of_propagators 1\n')
    input_file.write('\n')
    input_file.write('#Propagator {0}\n'.format(prop_no))                 
    input_file.write('mass {0}\n'.format(data['daughter prop']['mass']))
    input_file.write('naik_term_epsilon {0}\n'.format(naik(data['daughter prop']['mass'],data['daughter prop']['quark'])))
    input_file.write('error_for_propagator {0}\n'.format(data['daughter prop']['error']))
    input_file.write('rel_error_for_propagator {0}\n'.format(data['daughter prop']['rel_error']))
    if save == False:
        input_file.write('reload_serial_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist,spin_taste,data['lattice info']['tag']) )
    elif spin_taste != 'G5-G5':
        input_file.write('reload_serial_ksprop ./temp/{6}{0}.{1}_t{2}_wallprop_m{3}_tw{4}_st{5}\n'.format(data['lattice info']['ens'],cfg,t0,data['daughter prop']['mass'],twist,spin_taste,data['lattice info']['tag']) )
    else:
        input_file.write('reload_serial_ksprop ./props/{0}.{1}_wallprop_m{2}_th{3}_t{4}\n'.format(data['lattice info']['ens'],cfg,data['daughter prop']['mass'],twist,t0))
    input_file.write('forget_ksprop\n')
    input_file.write('\n')    
    return()

def no_sets_mesons_3pts(data):
    if data['parent prop']['multimass'] == True:
        set_no = len(data['parent prop']['spin_taste']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists'])  
    if data['parent prop']['multimass'] == False:
        set_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists'])  
        
    p_no = len(data['parent prop']['spin_taste'])*len(data['parent prop']['masses']) + len(data['daughter prop']['spin_taste'])*len(data['daughter prop']['twists']) 
        
    meson_no =  len(data['three points']['p J d'])*len(data['daughter prop']['twists'])*len(data['parent prop']['masses'])
    return(set_no,meson_no,p_no)


def make_mesons_3pt(data,filename,input_file,mass1,mass2,prop1,prop2,twist,t0,spin_taste,label):   
    input_file.write('\n')
    input_file.write('# ==================\n')
    input_file.write('# current-{5} masses {0} {1} {2} twist {3} {4}\n'.format(mass1,mass2,data['spectator prop']['mass'],twist,spin_taste,label))
    input_file.write('pair {0} {1}\n'.format(prop1,prop2))
    input_file.write('spectrum_request meson\n')
    input_file.write('save_corr_fnal {0}\n'.format(filename))
    input_file.write('r_offset 0 0 0 {0}\n'.format(t0))
    input_file.write('number_of_correlators 1\n')
    input_file.write('correlator ps p000 1 / {0} {1} 0 0 0 E E E\n'.format(data['lattice info']['wpnorm'],spin_taste))
    return()

def get_currents(data):
    currents = collections.OrderedDict()
    currents['p'] = []
    currents['d'] = []
    currents['J'] = []
    for i in range(len(data['three points']['p J d'])):
        currents['p'].append(data['three points']['p J d'][i].split()[0])
        currents['J'].append(data['three points']['p J d'][i].split()[1])
        currents['d'].append(data['three points']['p J d'][i].split()[2])       
    return(currents)


def main_3pts(argv):
    data = load_data()
    cfg = int(argv[0])
    input_file = open('./in/input-3pt/{0}milc_3pt_{1}.in'.format(data['lattice info']['tag'],cfg), 'w+')
    t0s = times(data,cfg)
    make_preamble_3pt(data,input_file,cfg)
    i = 0
    set_no,meson_no,prop_no = no_sets_mesons_3pts(data)
    currents = get_currents(data)
    for t0num,t0 in enumerate(t0s):
        Ts,actual_Ts = Times(data,cfg,t0)
        for j,T in enumerate(Ts):
            Props = collections.OrderedDict()
            Props['st'] = []
            Props['mass'] = []
            Props['twist'] = []
            Props['type'] = []
            linebreak('Source time {0} ({2} of {3}) , meson separation T = {1}'.format(t0,actual_Ts[j]-t0, t0num + 1,len(t0s)) ,input_file,80)
            make_gauge_field(data,cfg,input_file,i)
            i += 1 
            linebreak('Description of base sources',input_file,40)
            ###################### BASE SOURCES ################################
            input_file.write('number_of_base_sources {0}\n'.format(1+len(data['parent prop']['spin_taste'])))
            make_base_source('vec_prop',0,input_file,data,t0,cfg,'no_file_name_needed')
            for n,st in enumerate(data['parent prop']['spin_taste']):
                fname = './temp/{6}{0}.{1}_t{2}_extsrc_{3}_T{4}_m{5}'.format(data['lattice info']['ens'],cfg,t0,st,T,data['spectator prop']['mass'],data['lattice info']['tag'])
                make_base_source('vec_field',n+1,input_file,data,t0,cfg,fname)
            linebreak('Description of modified sources',input_file,40)
            ###################### MODIFIED SOURCES ############################
            input_file.write('number_of_modified_sources 0\n')
            linebreak('Description of propagators',input_file,40)
            ###################### PROPAGATORS #################################
            pr_num = 0
            set_num = 0
            input_file.write('number_of_sets {0}\n'.format(set_no))            
            #~~~~~~~~~~~~~~~~~~ DAUGHTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for st in remove_duplicates(data['daughter prop']['spin_taste']):
                for twist in data['daughter prop']['twists']:
                    if data['spectator prop']['same'] == True and st == 'G5-G5' and twist == data['spectator prop']['twist']:
                        if data['spectator prop']['save'] == True and data['spectator prop']['load'] == False:
                            saved = True
                        else:
                            saved = False
                        make_daughter_set_prop_3pt(Props,saved,input_file,data,twist,set_num,pr_num,0,st,cfg,t0)
                    else:
                        make_daughter_set_prop_3pt(Props,data['daughter prop']['save'],input_file,data,twist,set_num,pr_num,0,st,cfg,t0)
                    pr_num += 1 
                    set_num += 1
            #~~~~~~~~~~~~~~~~~~~~~ PARENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if data['parent prop']['multimass']== False:
                for st in remove_duplicates(data['parent prop']['spin_taste']):
                    for mass in data['parent prop']['masses']:
                        make_parent_set(Props,input_file,data,set_num,1+data['parent prop']['spin_taste'].index(st),st,1)
                        make_parent_prop(Props,input_file,data,mass,pr_num)
                        pr_num+=1
                        set_num+=1
            else:
                for st in remove_duplicates(data['parent prop']['spin_taste']):
                    make_parent_set(Props,input_file,data,set_num,1+data['parent prop']['spin_taste'].index(st),st,len(data['parent prop']['masses']))
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
            input_file.write('number_of_mesons {0}\n'.format(meson_no))
            for k in range(len(data['three points']['p J d'])):
                for l in range(len(Props['type'])): #l is daughter
                    for m in range(len(Props['type'])): # m is parent
                        if Props['type'][l] == 'daughter' and Props['st'][l] == currents['d'][k] and Props['type'][m] == 'parent' and Props['st'][m] == currents['p'][k]:
                            mass1 = Props['mass'][l]
                            mass2 = Props['mass'][m]
                            prop1 = int(l)
                            prop2 = int(m)
                            twist = Props['twist'][l]
                            spin_taste = currents['J'][k]
                            fname = './corrs/current-{0}_tw{1}/{9}current-{0}.{2}.{3}_t{4}_T{5}_m{6}_m{7}_m{8}_tw{1}'.format(data['three points']['label'][k],twist,data['lattice info']['ens'],cfg,t0,actual_Ts[j]-t0,mass1,mass2,data['spectator prop']['mass'],data['lattice info']['tag'])
                            make_mesons_3pt(data,fname,input_file,mass1,mass2,prop1,prop2,twist,t0,spin_taste,data['three points']['label'][k])
                            
            linebreak('Description of baryons',input_file,40)
            ################### BARYONS ##########################################
            input_file.write('number_of_baryons 0\n')
    input_file.close()            
    return()









def edit_extract(yamldata):    # edits extract/Extract.py
    if not os.path.exists('./extract'):
        os.makedirs('./extract')
    g = open('./extract/{0}Extract.py'.format(yamldata['lattice info']['tag']), 'w+')
    t0s = []
    Twists = []
    Masses = []
    Data = []
    twists = []
    masses = []
    for tnum in range(yamldata['lattice info']['nsrc']):
        t0s.append(tnum)
    if 'daughter prop' in yamldata:
        Data.append('{0}'.format(yamldata['daughter prop']['name']))
        for i,tw in enumerate(yamldata['daughter prop']['twists']):
            Twists.append(i)
            twists.append('{0}'.format(tw))
    if 'parent prop' in yamldata:
        Data.append('{0}'.format(yamldata['parent prop']['name']))
        for j,m in enumerate(yamldata['parent prop']['masses']):
            Masses.append(j)
            masses.append('{0}'.format(m))
    if yamldata['lattice info']['justtwopoints'] == False:
        for label in yamldata['three points']['label']:
           Data.append('{0}'.format(label))
    data = []
    if 'parent prop' in yamldata:
        for st in yamldata['parent prop']['spin_taste']:
            data.append('../corrs/{0}_{1}_tw{2}'.format(yamldata['parent prop']['name'],st,yamldata['parent prop']['twist']))
    if 'daughter prop' in yamldata:
        for st in yamldata['daughter prop']['spin_taste']:
            for twist in yamldata['daughter prop']['twists']:
                data.append('../corrs/{0}_{1}_tw{2}'.format(yamldata['daughter prop']['name'],st,twist))
    if yamldata['lattice info']['justtwopoints'] == False:
        for k in range(len(yamldata['three points']['p J d'])):
            for twist in yamldata['daughter prop']['twists']:
                data.append('../corrs/current-{0}_tw{1}'.format(yamldata['three points']['label'][k],twist))

    Ts = []
    for T in range(yamldata['three points']['nT']):
        Ts.append(yamldata['three points']['Tstart']+T*yamldata['three points']['dT'])
    g.write("import collections\n")
    g.write("import os.path\n")
    g.write("\n")
    g.write("# To edit\n")
    g.write("##########################################\n")
    g.write("cfgstart =\n")
    g.write("cfgend =\n")
    g.write("dconf =\n")
    g.write("skip = [] # any confs to not include\n")
    g.write("imaginary = [] # any imaginary corrs\n")
    g.write("t0s = {0}\n".format(t0s))
    g.write("Twists = {0}\n".format(Twists))
    g.write("Masses = {0}\n".format(Masses))
    g.write("Ts = {0}\n".format(Ts))
    g.write("Data = {0}\n".format(Data))
    g.write("negative = False\n")
    g.write("##########################################\n")
    g.write("\n")
    g.write("\n")
    g.write("\n")
    g.write("\n")
    g.write("twists = {0}\n".format(twists))
    g.write("masses = {0}\n".format(masses))
    g.write("data = {0}\n".format(data))
    g.write("lat ='{0}'\n".format(yamldata['lattice info']['ens']))
    g.write("nx = {0}\n".format(yamldata['lattice info']['nx']))
    g.write("ny = {0}\n".format(yamldata['lattice info']['ny']))
    g.write("nz = {0}\n".format(yamldata['lattice info']['nz']))
    g.write("nt = {0}\n".format(yamldata['lattice info']['nt']))
    g.write("norm = nx*ny*nz\n")
    g.write("mspec = {0}\n".format(yamldata['spectator prop']['mass']))
    g.write("mdaughter = {0}\n".format(yamldata['daughter prop']['mass']))
    g.write("nsrc = {0}\n".format(yamldata['lattice info']['nsrc']))
    g.write("dsrc = nt/nsrc\n")
    g.write("tag = '{0}'\n".format(yamldata['lattice info']['tag']))
    g.write("\n")
    g.write("labforgpl = Data[0]\n")
    g.write("for i,element in enumerate(Data):\n")
    g.write("    if i !=0:\n")
    g.write("        labforgpl ='{0}{1}'.format(labforgpl,element)\n")
    g.write("\n")
    g.write("gpl = '{0}_{1}cfgs_neg{2}.gpl'.format(labforgpl,int((cfgend-cfgstart)/dconf +1-len(skip)),negative)\n")
    g.write("print('Writing to {0}'.format(gpl))\n")
    g.write("\n")
    g.write("if os.path.exists(gpl):\n")
    g.write("    os.remove(gpl)\n")
    g.write("\n")
    g.write("def main(gpl):\n")
    g.write("    g = open(gpl,'w+')\n")
    g.write("    for cfg in range(cfgstart,cfgend+dconf,dconf):\n")
    g.write("        if cfg in skip:\n")
    g.write("            continue\n")
    g.write("        sources = []\n")
    g.write("        for n in range(len(t0s)):\n")
    g.write("            sources.append(int((19*(cfg/{0})) % dsrc + t0s[n]*dsrc))\n".format(yamldata['lattice info']['nsrcdivider']))
    if 'parent prop' in yamldata:
        g.write("        if '{0}' in Data:\n".format(yamldata['parent prop']['name']))
        g.write("            for element in data:\n")
        g.write("                if element.split('/')[2].split('_')[0] == '{0}':\n".format(yamldata['parent prop']['name']))
        g.write("                    for Mass in Masses:\n")
        g.write("                        mass = masses[Mass]\n")
        g.write("                        spin_taste = element.split('/')[2].split('_')[1]\n")
        g.write("                        filenames = []\n")
        g.write("                        datatag = '{0}_{1}_m{2}'.format(element.split('/')[2].split('_')[0],spin_taste,mass)\n")
        g.write("                        for source in sources:\n")
        g.write("                            filenames.append('{8}/{0}meson-{7}-{6}.{1}.{2}_t{3}_m{4}_m{5}'.format(tag,lat,cfg,source,mspec,mass,spin_taste,element.split('/')[2].split('_')[0],element))\n")
        g.write("                        if '{0}' in imaginary:\n".format(yamldata['parent prop']['name']))
        g.write("                            average_imaginary(filenames,datatag,g)\n")
        g.write("                        else:\n")
        g.write("                            average(filenames,datatag,g)\n")        
    if 'daughter prop' in yamldata:
        g.write("        if '{0}' in Data:\n".format(yamldata['daughter prop']['name']))
        g.write("            for element in data:\n")
        g.write("                if element.split('/')[2].split('_')[0] == '{0}':\n".format(yamldata['daughter prop']['name']))
        g.write("                    twist = element.split('/')[2].split('_tw')[1]\n")
        g.write("                    spin_taste = element.split('/')[2].split('_')[1]\n")
        g.write("                    filenames = []\n")
        g.write("                    datatag = '{0}_{1}_tw{2}'.format(element.split('/')[2].split('_')[0],spin_taste,twist)\n")
        g.write("                    for source in sources:\n")
        g.write("                        filenames.append('{8}/{0}meson-{7}-{6}.{1}.{2}_t{3}_m{4}_m{5}'.format(tag,lat,cfg,source,mspec,mdaughter,spin_taste,element.split('/')[2].split('_')[0],element))\n")
        g.write("                    if '{0}' in imaginary:\n".format(yamldata['daughter prop']['name']))
        g.write("                        average_imaginary(filenames,datatag,g)\n")
        g.write("                    else:\n")
        g.write("                        average(filenames,datatag,g)\n") 
    if yamldata['lattice info']['justtwopoints'] == False:
        for threeptlabel in yamldata['three points']['label']:
            g.write("        if '{0}' in Data:\n".format(threeptlabel))
            g.write("            for element in data:\n")
            g.write("                if element.split('/')[2].split('-')[1].split('_')[0] == '{0}':\n".format(threeptlabel))
            g.write("                    twist = element.split('/')[2].split('_tw')[1]\n")
            g.write("                    for Mass in Masses:\n")
            g.write("                        mass = masses[Mass]\n")
            g.write("                        for T in Ts:\n")
            g.write("                            filenames = []\n")
            g.write("                            datatag = '{0}_T{1}_m{2}_m{3}_m{4}_tw{5}'.format(element.split('/')[2].split('-')[1].split('_')[0],T,mdaughter,mass,mspec,twist)\n")
            g.write("                            for source in sources:\n")
            g.write("                                filenames.append('{0}/{10}current-{1}.{3}.{4}_t{5}_T{6}_m{7}_m{8}_m{9}_tw{2}'.format(element,element.split('/')[2].split('-')[1].split('_')[0],twist,lat,cfg,source,T,mdaughter,mass,mspec,tag))\n")
            g.write("                            if '{0}' in imaginary:\n".format(threeptlabel))
            g.write("                                average_imaginary(filenames,datatag,g)\n")
            g.write("                            else:\n")
            g.write("                                average(filenames,datatag,g)\n") 
    g.write("        print('Extracted conf', cfg)\n")
    g.write("    g.close()\n")
    g.write("    return()\n")
    g.write("\n")
    g.write("\n")
    g.write("def test_zeros(filename):\n")
    g.write("    if os.stat(filename).st_size == 0:\n")
    g.write("        print(filename, 'empty')\n")
    g.write("        return(0,0)\n")
    g.write("    f = open(filename,'r')\n")
    g.write("    lastDots = 0\n")
    g.write("    lines = f.readlines()\n")
    g.write("    for i, line in enumerate(lines):\n")
    g.write("        if line == '...\\n':\n")
    g.write("            lastDots = i\n")
    g.write("    if lines[lastDots + 5].split()[1] == '0.000000e+00':\n")
    g.write("        print(filename, 'contains zeros')\n")
    g.write("    f.close()\n")
    g.write("    return(lines,lastDots)\n")
    g.write("\n")
    g.write("\n")
    g.write("def average(filenames,datatag,g):\n")
    g.write("    result = collections.OrderedDict()   \n")
    g.write("    for filename in filenames:                          #filenames is list of files to be averaged\n")
    g.write("        result[filename] = []\n")
    g.write("        if os.path.exists(filename):\n")
    g.write("            lines,lastdots = test_zeros(filename)\n")
    g.write("            if lines == 0:\n")
    g.write("                continue\n")
    g.write("        else:\n")
    g.write("            print(filename, 'Does not exist')\n")
    g.write("            continue\n")
    g.write("        for j in range(lastdots+1,lastdots+1+nt):\n")
    g.write("            result[filename].append(lines[j].split()[1])\n")
    g.write("    g.write('{0}    '.format(datatag))\n")
    g.write("    for i in range(nt):\n")
    g.write("        average = 0\n")
    g.write("        for filename in filenames:\n")
    g.write("            if os.path.exists(filename):\n")
    g.write("                average += float(result[filename][i])/(norm*len(filenames))\n")
    g.write("        if negative == True:\n")
    g.write("            average = -average\n")
    g.write("        g.write('{0:.11g}    '.format(average))\n")
    g.write("    g.write('\\n')\n")
    g.write("    return()\n")
    g.write("\n")
    g.write("\n")
    g.write("def average_imaginary(filenames,datatag,g):\n")
    g.write("    result = collections.OrderedDict()   \n")
    g.write("    for filename in filenames:                          #filenames is list of files to be averaged\n")
    g.write("        result[filename] = []\n")
    g.write("        if os.path.exists(filename):\n")
    g.write("            lines,lastdots = test_zeros(filename)\n")
    g.write("            if lines == 0:\n")
    g.write("                continue\n")
    g.write("        else:\n")
    g.write("            print(filename, 'Does not exist')\n")
    g.write("            continue\n")
    g.write("        for j in range(lastdots+1,lastdots+1+nt):\n")
    g.write("            result[filename].append(lines[j].split()[2])\n")
    g.write("    g.write('{0}    '.format(datatag))\n")
    g.write("    for i in range(nt):\n")
    g.write("        average = 0\n")
    g.write("        for filename in filenames:\n")
    g.write("            if os.path.exists(filename):\n")
    g.write("                average += float(result[filename][i])/(norm*len(filenames))\n")
    g.write("        if negative == False:\n")
    g.write("            average = -average\n")
    g.write("        g.write('{0:.11g}    '.format(average))\n")
    g.write("    g.write('\\n')\n")
    g.write("    return()\n")
    g.write("\n")
    g.write("\n")
    g.write("main(gpl)\n")
    g.write("\n")
    g.close()    
    return()



edit_extract(data)


main_2pts(sys.argv[1:])
if data['lattice info']['justtwopoints'] == False:
    main_extsrc(sys.argv[1:])
    main_3pts(sys.argv[1:])
