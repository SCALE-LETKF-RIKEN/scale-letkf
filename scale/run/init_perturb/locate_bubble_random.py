
import sys

def read_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--ensemble_size', type=int, required=True)
    parser.add_argument('-cx', '--center_of_bubble_x', type=float, required=True)
    parser.add_argument('-cy', '--center_of_bubble_y', type=float, required=True)
    parser.add_argument('-stdx', '--loc_std_x', type=float, required=True)
    parser.add_argument('-stdy', '--loc_std_y', type=float, required=True)
    parser.add_argument('-nml', '--namelist_basename', type=str, required=True)
    parser.add_argument('-npz', '--npz_directory', type=str, required=True)
    return parser.parse_args()

def calc_random_loc(ensemble_size, center_x, center_y, std_x, std_y, npz_dir):
    npz_file = '{0:}/loc_ens{1:05}_cx{2:.1f}_cy{3:.1f}_stdx{4:.1f}_stdy{5:.1f}.npz'.format(
        npz_dir, ensemble_size, center_x, center_y, std_x, std_y)
    try:
        loc_x, loc_y = read_loc_from_npz(npz_file)
        return loc_x, loc_y
    except:
        print('calculate random location')
        
    import numpy as np
    loc_x = np.random.normal(center_x, std_x, ensemble_size)
    loc_y = np.random.normal(center_y, std_y, ensemble_size)
    
    # Set the mean location
    loc_x = loc_x - np.mean(loc_x) + center_x
    loc_y = loc_y - np.mean(loc_y) + center_y 
    
    # Normalize the standard deviation
    loc_x = loc_x / np.std(loc_x) * std_x
    loc_y = loc_y / np.std(loc_y) * std_y
    
    import os
    os.makedirs(npz_dir, exist_ok=True)
    np.savez(npz_file, loc_x=loc_x, loc_y=loc_y, center_x=center_x, center_y=center_y, std_x=std_x, std_y=std_y)
    return loc_x, loc_y

def read_loc_from_npz(npz_file):
    import numpy as np
    data = np.load(npz_file)
    loc_x = data['loc_x']
    loc_y = data['loc_y']
    return loc_x, loc_y

def replace_bbl_loc_in_namelist(loc_x=None, loc_y=None, namelist=None):

    with open(namelist, 'r') as f:
        lines = f.readlines()

    with open(namelist, 'w') as f:
        for line in lines:
            if '!--BBL_CX--' in line:
                f.write(' BBL_CX = {0:.2f}d0,\n'.format(loc_x))
            elif '!--BBL_CY--' in line:
                f.write(' BBL_CY = {0:.2f}d0,\n'.format(loc_y))
            else:
                f.write(line)
    

def main():
    args = read_args()
    loc_x, loc_y = calc_random_loc(args.ensemble_size, args.center_of_bubble_x, args.center_of_bubble_y, 
                                   args.loc_std_x, args.loc_std_y, args.npz_directory)

    namelist_basename = args.namelist_basename
    
    for i in range(args.ensemble_size):
        namelist = namelist_basename.replace('fxxxx','f{0:04}'.format(i+1))
        replace_bbl_loc_in_namelist(loc_x[i], loc_y[i], namelist=namelist)
            
if __name__ == '__main__':
    main()