#! /usr/bin/env python
"""Convert ASTRA output file into an INPUT_profiles file for GYRO.

   Usage to call as a script:
     astra2gyro.py [options] path
     [options]: optional aguments:
     -o <output file>, --output <output file>: write output to <output file>
     -m <map file>, --map <map file>: use <map file> to map ASTRA data to GYRO INTPUT_profiles headers
                                      uses 'astra_map_gyro' in current directory by default
     -n <note text>, --note <note text>: include <note text> in the INPUT_profiles header
   Example usage:"
        astra2gyro.py -m ../base/astra_map_gyro -o INPUT_profiles -n 'this is a test' astra_data_file.dat

   Alternatively, can be imported as a module during an interactive python session for greater flexibility.
     For example:
       >>>import astra2gyro as ag
       >>>args=ag.defaultArguments('astra_data_file.dat')
       >>>out=ag.main(args)
       >>>out['gyro']['q']=map(lambda r:out['gyro']['q'][0]*
                               (r/out['gyro']['rmin'][-1])**2.0,out['gyro']['rmin'])
       >>>writeGyroInputProfiles(out['gyro'],out['header']+'set to a parabolic q profile\n','INPUT_profiles')

   The map file links data names in the ASTRA data file to the corresponding names for GYRO INPUT_profiles.
   See the included sample map file: 'astra_map_gyro'.
       
   See the docstrings for each provided function for more information.
"""
"""@package astra2gyro
   @author Erik Granstedt
   @date 2009-06-30
"""

# --- globals ---
astra_delimiter = ' '           # or '\t'
map_default='astra_map_gyro'
output_extension_default='.gyro'
neg_sign='-'
file_header="""#
# ** This file was generated by astra2gyro.py **
#
"""

INPUT_profiles_scalar=['BT_EXP','N_EXP','ARHO_EXP'] # scalars come first
INPUT_profiles_vector=[['rho','rmin','rmaj','q','kappa'],     # vector data block 1
                       ['delta','te','ne','z_eff','er'],
                       ['flow_mom','pow_e','pow_i','pow_ei','zeta'],
                       ['flow_beam','flow_wall','zmag'],
                       ['ni_1','ni_2','ni_3','ni_4','ni_5'],
                       ['ti_1','ti_2','ti_3','ti_4','ti_5'],
                       ['vtor_1','vtor_2','vtor_3','vtor_4','vtor_5'],
                       ['vpol_1','vpol_2','vpol_3','vpol_4','vpol_5']] # vector data block 8
INPUT_profiles_delimiter='   '
INPUT_profiles_defaults={'BT_EXP':1.0,'N_EXP':0,'ARHO_EXP':1.0, # scalars
                         'rho':0.0,'rmin':0.0,'rmaj':0.0,'q':0.0,'kappa':1.0,     # vector data block 1
                         'delta':0.0,'te':0.0,'ne':0.0,'z_eff':1.0,'er':0.0,
                         'flow_mom':0.0,'pow_e':0.0,'pow_i':0.0,'pow_ei':0.0,'zeta':0.0,
                         'flow_beam':0.0,'flow_wall':0.0,'zmag':0.0,
                         'ni_1':0.0,'ni_2':0.0,'ni_3':0.0,'ni_4':0.0,'ni_5':0.0,
                         'ti_1':0.0,'ti_2':0.0,'ti_3':0.0,'ti_4':0.0,'ti_5':0.0,
                         'vtor_1':0.0,'vtor_2':0.0,'vtor_3':0.0,'vtor_4':0.0,'vtor_5':0.0,
                         'vpol_1':0.0,'vpol_2':0.0,'vpol_3':0.0,'vpol_4':0.0,'vpol_5':0.0} # vector data block 8


# --- modules to import
from re import search
from os.path import abspath

# --- function definitions
def has_letter(s):
    """Tests if the string contains any non-numerical characters."""
    if search(r'[A-df-z!@\$%\^&\*\(\)\[\]\{\}_=<>\"\'\?\\/]+',s): # note 'e' used in scientific notation ex: 1.005e-02
        return True
    else:
        return False

def processSingleValues(scalars):
    """Process description lines from beginning and end of ASTRA output.

    Returns a dictionary with variable names and their values, as well
    as a list of the remaining description strings.
    """
    vals=map(astraAssignedToFloat,scalars)
    desc=map(lambda x:x[0],filter(lambda x:len(x)==1,vals))
    vals=filter(lambda x:len(x)==2,vals)
    data={}
    for x in vals:
        data[x[0]]=x[1]
    return [data,desc]

def astraAssignedToFloat(x):
    """Convert a string of the form "string=value" into [string,value] pairs."""
    elems=x.split('=')
    if len(elems)==1:
        if elems[0]:
            try:
                y=astra_val_to_float(elems[0])
                return [y]
            except ValueError:
                return [elems[0]]
        else:
            return []
    if len(elems)==2:
        if elems[1]:
            try:
                y=astra_val_to_float(elems[1])
                return [elems[0],y]
            except ValueError:
                return [x]
        else:
            return [elems[0]]
    if len(elems)>2:
        return []

def astra_val_to_float(x):
    """Convert ASTRA strings like "967-6" to floats like 9.67e-6."""
    try:
        y=float(x)
        return y
    except ValueError:
        neg=False
        if x[0]=='-':           # test if negative
            x=x[1:]
            neg=True
        if x[0]=='+':           # test if explicit plus sign
            x=x[1:]
        if '+' in x:            # test if positive power of 10
            elems=x.split('+')
            if len(elems)==2:
                y=float(elems[0][0]+'.'+elems[0][1:])*10.**float(elems[1])
            else:
                raise ValueError, x
        elif '-' in x:          # test if neg. power of 10
            elems=x.split('-')
            if len(elems)==2:
                y=float(elems[0][0]+'.'+elems[0][1:])*10.**(-float(elems[1]))
            else:
                raise ValueError, x
        else:
            raise ValueError, x
        return -y if neg else y

def fix_spurious_neg_sign(string_arr):
    """Fix negative signs that have been separated from their corresponding numbers"""
    for i in range(string_arr.count(neg_sign)):
        idx=string_arr.index(neg_sign)
        string_arr.pop(idx)
        string_arr[idx]=neg_sign+string_arr[idx]

def astra_cols_to_list(flines):
    """Returns a list of rows as floating point values."""
    #flines=map(lambda line:line[2:-1],flines)  # lines always have 3 spaces at beginning, 1 at end
    #flines=map(lambda line:map(lambda i:line[1+i*7:(i+1)*7],range(len(line)/7)),flines)# each column is always 6 characters long
    
    flines=map(lambda line:line.split(),flines)            # split the data at separators
    flines=remove_whitespace_blanks(flines)                # remove space and blanks
    map(lambda line:fix_spurious_neg_sign(line),flines) # fix negative signs that have been separated from their numbers
    # convert values like "967-6" to 9.67e-6
    data=map(lambda line: map(astra_val_to_float,line),flines) # convert strings to list of floats
    return data

def remove_whitespace_blanks(list_of_strings):
    """Remove whitespace before and after each string in a list of strings; also remove empty strings."""
    x=map(lambda z:map(lambda y:y.strip(),z),list_of_strings) # remove whitespace
    x=map(lambda z:filter(lambda i: i != '',z),x) # remove empties
    return x

def readAstraFile(fname):
    """Read in an ASTRA output file; return as a dictionary with name:value[s]."""
    # idea: read in file
    # lines with an equal sign split into separate lines at the tabs
    # identify the lines that are headers
    # read in the data lines into their appropriate names
    f=open(fname,'r')
    flines=f.readlines()
    f.close()
    # process scalar data
    single_value_lines=filter(lambda x:True if '=' in x else False,flines)    # lines with equal signs
    single_values=sum(map(lambda x:x.split(),single_value_lines),[]) # flatten list
    single_values=map(lambda x:x.strip(),single_values) # remove whitespace
    single_values=filter(lambda i: i != '', single_values) # remove empties
    single_value_line_nums=map(lambda x:flines.index(x),single_value_lines)
    [scalars,description]=processSingleValues(single_values)     # split scalar data from text
    # header lines
    headers=filter(has_letter,flines[single_value_line_nums[0]+1:single_value_line_nums[1]])
    header_line_nums=map(lambda x:flines.index(x),headers)
    headers=map(lambda x:x.split(),headers) # split at delimiters
    cols_per_header_line=map(lambda x:len(x),headers)
    headers=sum(headers,[]) # flatten list
    headers=map(lambda x:x.strip(),headers) # remove whitespace
    headers=filter(lambda i: i != '', headers) # remove empties
    # lines with only numeric data
    data_arr_list=[]
    for i in range(len(header_line_nums)-1):
        theseflines=flines[header_line_nums[i]+1:header_line_nums[i+1]]
        data_arr_list.append(astra_cols_to_list(theseflines))
    theseflines=flines[header_line_nums[-1]+1:single_value_line_nums[1]]
    data_arr_list.append(astra_cols_to_list(theseflines))
    # create a dictionary matching the header name to the values
    data={}
    i=0
    for block in data_arr_list:
        for colj in range(len(block[0])):
            vals=map(lambda line:line[colj],block)
            data[headers[i]]=vals if len(vals) > 1 else vals[0]
            i=i+1
    return {'description':description,'vector':data,'scalar':scalars}#,data_arr_list,[single_value_line_nums,header_line_nums,cols_per_header_line]]

def processMapValue(key,value,astra_dict,gyro_dict,nradial):
    """Map a key in the map file to its astra data."""
    try:
        # scalar data
        if key in reduce(lambda x,y:x+y,INPUT_profiles_scalar):
            if value.has_key('override'):        # override to a given value
                gyro_dict[key]=value['override']
            elif value.has_key('astra'):# referencing an astra dictionary value
                if value.has_key('index'): # reference a particular index of an astra vector
                    gyro_dict[key]=astra_dict['vector'][value['astra']][value['index']]
                else:               # reference an astra scalar
                    gyro_dict[key]=astra_dict['scalar'][value['astra']]
            elif not value: # if empty dictionary choose default
                gyro_dict[key]=nradial if key=='N_EXP' else INPUT_profiles_defaults[key]
            else:               # dictionary has wrong form
                raise ValueError
        # vector data
        elif key in reduce(lambda x,y:x+y,INPUT_profiles_vector):
            if value.has_key('override'):        # override to a given value
                if len(value['override'])==nradial:
                    gyro_dict[key]=value['override']
                elif len(value['override'])==1:
                    gyro_dict[key]=[value['override']]*nradial
            elif value.has_key('astra'):# referencing an astra dictionary value
                if value.has_key('index'): # reference a particular index of an astra vector
                    gyro_dict[key]=astra_dict['vector'][value['astra']][value['index']]*nradial
                else:               # reference an entire astra vector
                    gyro_dict[key]=astra_dict['vector'][value['astra']]
            elif not value: # if empty dictionary choose default
                gyro_dict[key]=[INPUT_profiles_defaults[key]]*nradial
            else:               # dictionary has wrong form
                raise ValueError
        else:                   # key is not an INPUT_profiles scalar or vector
            raise ValueError
    except ValueError:
        print 'gyroFromAstra: key:%s in map file is not valid'%key
        raise ValueError    # re-throw exception and exit

def gyroFromAstra(astra_dict,mapfile):
    """Use file 'map' to map a dictionary of astra profile values into GYRO profile values."""
    map_elems={}
    execfile(mapfile,{},map_elems)
    gyro_dict={}
    map_elems.pop('__doc__')   # discard this key, comes from running with 'execfile'
    if map_elems.has_key('nonmap'):
        nonmap=map_elems.pop('nonmap') # get the list of keys that don't map to astra data
        for item in nonmap:
            # get the lists which give the order of elements in GYRO INPUT_profiles file
            # these keys do not map to astra data
            gyro_dict[item]=map_elems.pop(item)
    # determine number of radial gridpoints
    if map_elems['N_EXP']:
        processMapValue('N_EXP',map_elems['N_EXP'],astra_dict,gyro_dict,0)
        nradial=gyro_dict['N_EXP']
    elif map_elems['rho']:
        processMapValue('rho',map_elems['rho'],astra_dict,gyro_dict,0)
        nradial=len(gyro_dict['rho'])
    elif map_elems['rmin']:
        processMapValue('rmin',map_elems['rmin'],astra_dict,gyro_dict,0)
        nradial=len(gyro_dict['rmin'])
    else:
        print "gyroFromAstra: number of radial points not specified through either, 'N_EXP', 'rho', or 'rmin'."
        raise ValueError
    # map other values from map file to the gyro INPUT_profiles dictionary
    for key in map_elems:
        processMapValue(key,map_elems[key],astra_dict,gyro_dict,nradial)
 
    return gyro_dict

def makeHeader(fname,mapfile,note=""):
    """Return header string for the GYRO INPUT_profiles file."""
    header=file_header+"#   ASTRA FILE: %s\n#   MAP FILE: %s\n"%(abspath(fname),abspath(mapfile))
    header=header+"#   note: %s\n"%note if note else header
    return header

def transpose(data):
    """Return transpose of 2D list."""
    return map(lambda i:map(lambda row:row[i],data),range(len(data[0])))

def writeGyroInputProfiles(gyro_dict,header,fname):
    """Write GYRO INPUT_profiles file to 'fname' from a dictionary of values, 'gyro_dict'."""
    f=open(fname,'w')
    f.write(header+'#\n')
    data_block_length=len(INPUT_profiles_vector[0])
    # write scalar data
    for item in INPUT_profiles_scalar:
        f.write(item+'='+str(gyro_dict[item])+'\n')
    # write vector data blocks
    for block in INPUT_profiles_vector:
        data=[]
        for item in block:
            data+=[gyro_dict[item]]
            column_headers=reduce(lambda x,y:'%12s'%x+INPUT_profiles_delimiter+'%13s'%y,block)
        if len(data)<data_block_length:
            data+=[[0]*len(data[0])]*(data_block_length-len(data))
        data=transpose(data)
        f.write('#\n#'+column_headers+'\n')
        data0_str=reduce(lambda str_elems,elem:str_elems+INPUT_profiles_delimiter+'%10.7e'%elem,['%10.7e'%data[0][0]]+data[0][1:])+'\n'
        data_str=reduce(lambda str_lines,line: str_lines+reduce(lambda str_elems,elem:str_elems+INPUT_profiles_delimiter+'%10.7e'%elem,['%10.7e'%line[0]]+line[1:])+'\n',[data0_str]+data[1:])
        f.write(data_str)
        #savetxt(f,data,fmt='%10.7e',delimiter=INPUT_profiles_delimiter)
    f.close()
    return 0

def main(arguments):
    """Produce a GYRO INPUT_profiles file from ASTRA output."""
    # read in datafiles
    fname = arguments['fname']
    mapfile = arguments['map'] if arguments['map'] else map_default
    output = arguments['output'] if arguments['output'] else fname+output_extension_default
    try:                        # read in the astra data file
        astra_dict=readAstraFile(fname)
    except IOError:
        print "astra2gyro: File '%s' not found."%fname
        return []
    gyro_dict=gyroFromAstra(astra_dict,mapfile) # map the astra data to gyro INPUT_profiles data
    header=makeHeader(fname,mapfile,note=arguments['note'])
    writeGyroInputProfiles(gyro_dict,header,output)
    return {'astra':astra_dict,'gyro':gyro_dict,'header':header}

def defaultArguments(fname,output='',mapfile='',note=''):
    """Create dictionary with the default aguments."""
    return {'fname':fname,'output':output,'map':mapfile,'note':note}

# --- script ---
if __name__ == "__main__":  #executing as script
    import sys, getopt
    # examine arguments; by default do not save figure, and do not run in batch mode
    arguments={'output':'','map':'','note':''}
    opts, extraparams = getopt.getopt(sys.argv[1:],'o:m:n:',['output=','map=','note='])
    for o,p in opts:
        if o in ['-o','--output']:
            arguments['output']=p
        if o in ['-m','--map']:
            arguments['map']=p
        if o in ['-n','--note']:
            arguments['note']=p

    if len(extraparams) != 1:
        print "Usage:"
        print "astra2gyro.py [options] path"
        print "[options]: optional aguments:"
        print "  -o <output file>, --output <output file>: write output to <output file>"
        print "  -m <map file>, --map <map file>: use <map file> to map Astra data to GYRO INTPUT_profiles headers"
        print "                                   uses 'astra_map_gyro' in current directory by default"
        print "  -n <note text>, --note <note text>: include <note text> in the INPUT_profiles header"
        print "example usage:"
        print "astra2gyro -m astra_map_gyro -o INPUT_profiles -n 'this is a test' LTXn.zbltx1.1"
        sys.exit(2)

    arguments['fname']=extraparams[0]
    a=main(arguments)
#    sys.exit(0 if main(arguments) else 1)

