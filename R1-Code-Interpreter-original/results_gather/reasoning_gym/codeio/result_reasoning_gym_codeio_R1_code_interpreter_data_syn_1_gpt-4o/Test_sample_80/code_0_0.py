import numpy as np

def compute_refractive_index(f, res, nm, method):
    f = np.array(f)
    
    if method == "ODT":
        km = (2 * np.pi * nm) / res
        ri = nm * np.sqrt(f / km**2 + 1)
        negrootcoord = np.where(ri.real < 0)
        ri[negrootcoord] *= -1
    elif method == "OPT":
        ri = nm + f / (2 * np.pi) * res
    else:
        raise ValueError("Invalid method. Choose either 'ODT' or 'OPT'.")
    
    return ri.tolist()

# Given input
input_data = {
    'f': [5.840031058612119, 3.0586858609653085, 4.703237071129676, 4.267385316960167, 6.669199490618037, 3.688666637146013, 5.31025989057464, 0.3687838253547492, 5.9944621915622145, 2.4384947256322835, 4.614861189104287, 2.700731637926814, 2.032500337444157, 7.362909325716278, 6.768552397536424],
    'res': 3.594280244804947,
    'nm': 1.7465588633140519,
    'method': 'OPT'
}

# Compute the result
result = compute_refractive_index(input_data['f'], input_data['res'], input_data['nm'], input_data['method'])
print(result)