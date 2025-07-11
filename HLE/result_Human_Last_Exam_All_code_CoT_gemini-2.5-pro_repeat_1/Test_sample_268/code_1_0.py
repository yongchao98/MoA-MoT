import numpy as np

def count_non_compact_roots():
    """
    Calculates the number of non-compact positive roots for a real form of C8
    specified by a Vogan diagram.
    """
    # Rank of the Lie algebra
    n = 8

    # Vogan diagram: W -- B -- W -- B -- B -- W -- B == B
    # Simple roots are alpha_1, ..., alpha_8.
    # Black (non-compact) simple roots are {alpha_2, alpha_4, alpha_5, alpha_7, alpha_8}.
    # We use 0-based indexing for arrays, so the indices are {1, 3, 4, 6, 7}.
    black_indices = {1, 3, 4, 6, 7}

    def get_coeffs_in_simple_basis(root_in_e_basis):
        """
        Converts a root from the standard {e_i} basis to the simple root {alpha_i} basis.
        
        The transformation is based on:
        e_i = sum_{k=i to n-1} alpha_k + 0.5 * alpha_n  (for i < n)
        e_n = 0.5 * alpha_n
        
        A root R = sum(r_i * e_i) has coefficients c_j in the alpha basis given by:
        c_j = sum_{i=1 to j} r_i  (for j < n)
        c_n = 0.5 * sum_{i=1 to n} r_i
        """
        coeffs = np.zeros(n)
        # Calculate coefficients for alpha_1 to alpha_{n-1}
        for j in range(n - 1):
            # c_j = sum of r_i for i <= j
            coeffs[j] = np.sum(root_in_e_basis[:j+1])
        
        # Calculate coefficient for alpha_n
        coeffs[n - 1] = 0.5 * np.sum(root_in_e_basis)
        
        return coeffs

    # Counters for the two types of positive roots
    non_compact_type1_count = 0  # For roots of type e_i - e_j
    non_compact_type2_count = 0  # For roots of type e_i + e_j

    # 1. Iterate through positive roots of type e_i - e_j (for 1 <= i < j <= n)
    for i in range(n):
        for j in range(i + 1, n):
            root_e = np.zeros(n)
            root_e[i] = 1
            root_e[j] = -1
            
            coeffs = get_coeffs_in_simple_basis(root_e)
            
            black_coeff_sum = sum(coeffs[k] for k in black_indices)
            
            # Check if the sum is an odd integer (within a small tolerance)
            if abs(black_coeff_sum - round(black_coeff_sum)) < 1e-9 and round(black_coeff_sum) % 2 == 1:
                non_compact_type1_count += 1

    # 2. Iterate through positive roots of type e_i + e_j (for 1 <= i <= j <= n)
    for i in range(n):
        for j in range(i, n):
            root_e = np.zeros(n)
            root_e[i] = 1
            root_e[j] += 1 # Handles i=j case (2*e_i)
            
            coeffs = get_coeffs_in_simple_basis(root_e)
            
            black_coeff_sum = sum(coeffs[k] for k in black_indices)

            # Check if the sum is an odd integer
            if abs(black_coeff_sum - round(black_coeff_sum)) < 1e-9 and round(black_coeff_sum) % 2 == 1:
                non_compact_type2_count += 1

    # Print the breakdown of the calculation as an "equation"
    total_non_compact = non_compact_type1_count + non_compact_type2_count
    print(f"Number of non-compact positive roots of type e_i-e_j: {non_compact_type1_count}")
    print(f"Number of non-compact positive roots of type e_i+e_j: {non_compact_type2_count}")
    print(f"Total number of non-compact positive roots = {non_compact_type1_count} + {non_compact_type2_count} = {total_non_compact}")

    return total_non_compact

# Execute the function to get the final answer
final_answer = count_non_compact_roots()
# The final answer is wrapped below as requested
# print(f"\n<<< {final_answer} >>>")