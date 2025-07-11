import math
import numpy as np

def solve_vest_b(v, T_matrices, S, k):
    """
    Solves the VEST problem for the specific case in question (b).
    
    Args:
        v (list[float]): The vector v.
        T_matrices (list[list[list[float]]]): The list of matrices T_i.
        S (list[list[float]]): The matrix S.
        k (int): The parameter k.
    
    Returns:
        float: The final computed value.
    """
    n = len(v)
    m = len(T_matrices)
    
    # C_counts[c] will store the number of T_i matrices that equal E_cc
    C_counts = [0] * n
    
    for T in T_matrices:
        # Check if T is a valid E_cc matrix
        non_zero_pos = -1
        non_zero_count = 0
        is_valid = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    if T[i][j] == 1:
                        non_zero_pos = i
                        non_zero_count += 1
                    elif T[i][j] != 0:
                        is_valid = False
                        break
                elif T[i][j] != 0:
                    is_valid = False
                    break
            if not is_valid:
                break
        
        if is_valid and non_zero_count == 1:
            C_counts[non_zero_pos] += 1
            
    total_value = 0
    equation_parts = []
    
    for c in range(n):
        try:
            # Calculate binomial coefficient comb(|C_c|, k)
            num_choices = math.comb(C_counts[c], k)
        except ValueError:
            num_choices = 0
            
        # v_c^2 * S_cc
        value_per_choice = (v[c]**2) * S[c][c]
        
        term_value = num_choices * value_per_choice
        total_value += term_value
        
        # Storing parts for the final equation output
        part_str = f"{num_choices} * ({v[c]}^2) * {S[c][c]}"
        equation_parts.append(part_str)

    # Print the full calculation as an equation
    equation_str = " + ".join(equation_parts)
    print("The calculation is based on the formula: sum_c (comb(|C_c|, k) * v_c^2 * S_cc)")
    print("Final Equation:")
    print(f"{equation_str} = {total_value}")

    return total_value

# Example Usage:
n_dim = 4
k_param = 3
v_vec = [1, 2, 3, 4]
S_mat = np.array([[0, 1, 2, 3], [1, 1, 2, 3], [2, 2, 4, 6], [3, 3, 6, 9]]).tolist()

# T_matrices are Z_2 diagonal matrices with at most one non-zero entry
T_mats = [
    # 4 matrices are E_00
    [[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
    [[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
    [[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
    [[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
    # 2 matrices are E_11
    [[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]],
    [[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]],
    # 3 matrices are E_22
    [[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]],
    [[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]],
    [[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]],
    # 0 matrices are E_33
    # 1 zero matrix
    [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
]

solve_vest_b(v_vec, T_mats, S_mat, k_param)