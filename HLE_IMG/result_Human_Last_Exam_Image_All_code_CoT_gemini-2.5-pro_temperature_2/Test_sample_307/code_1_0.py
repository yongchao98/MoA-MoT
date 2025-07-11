import sympy
from sympy.physics.wigner import wigner_3j

def calculate_max_abs_3j(j1, j2, j3):
    """
    Calculates the maximum absolute value of the Wigner 3j symbol for a given
    (j1, j2, j3) over all valid m1 and m2 projections.
    This corresponds to the infinity-norm.
    """
    max_abs_val = 0.0
    
    # Iterate through all allowed m quantum numbers.
    # m1 ranges from -j1 to j1, m2 ranges from -j2 to j2.
    for m1_int in range(-int(j1), int(j1) + 1):
        for m2_int in range(-int(j2), int(j2) + 1):
            # The third m-number is fixed by the selection rule m1 + m2 + m3 = 0.
            m3_int = -m1_int - m2_int
            
            # The symbol is non-zero only if |m3| <= j3.
            if abs(m3_int) > j3:
                continue

            # Calculate the Wigner 3-j symbol using sympy.
            # It's good practice to use sympy's Integer types.
            val = wigner_3j(j1, j2, j3, m1_int, m2_int, m3_int)
            
            # The result from sympy can be a symbolic expression. Evaluate it to a float.
            abs_val_f = abs(val.evalf())
            
            if abs_val_f > max_abs_val:
                max_abs_val = abs_val_f
                
    return max_abs_val

def find_wigner_norm_ratio():
    """
    Identifies the parameters for the nine Wigner 3-j symbols from the image,
    calculates their infinity-norms, and finds the ratio of the max to min norm.
    """
    # Based on visual analysis and known visualizations, the nine symbols correspond to
    # configurations with j3=10 and j1, j2 chosen from {8, 9, 10}.
    # We create a list of these 9 (j1, j2, j3) configurations.
    j_configurations = []
    j3_fixed = 10
    for j1 in [8, 9, 10]:
        for j2 in [8, 9, 10]:
            # The triangle inequality |j1-j2| <= j3 <= j1+j2 must be satisfied.
            if abs(j1 - j2) <= j3_fixed <= j1 + j2:
                 j_configurations.append((j1, j2, j3_fixed))

    # Calculate the infinity-norm for each configuration.
    norms = []
    for params in j_configurations:
        j1, j2, j3 = params
        norm = calculate_max_abs_3j(j1, j2, j3)
        norms.append(norm)

    # Find the maximum and minimum norms from the list.
    if not norms:
        print("No valid Wigner 3-j symbol configurations found.")
        return

    max_norm = max(norms)
    min_norm = min(norms)

    # Calculate the ratio.
    ratio = max_norm / min_norm
    
    print("This script calculates the ratio of the maximum to the minimum infinity-norm among a set of nine Wigner 3-j symbols.")
    print("Based on analysis of the provided image, the symbols are for j3=10, with j1 and j2 varying in {8, 9, 10}.\n")
    print(f"Maximum infinity-norm found: {max_norm}")
    print(f"Minimum infinity-norm found: {min_norm}\n")
    print("The final equation for the ratio is:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    find_wigner_norm_ratio()
