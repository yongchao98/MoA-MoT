import numpy as np

def solve_vest_b(v, T_matrices, S, k):
    """
    Solves the VEST problem for the special case in question (b).
    
    Args:
        v (np.ndarray): The rational vector v.
        T_matrices (list of np.ndarray): The list of matrices T_i.
        S (np.ndarray): The matrix S.
        k (int): The parameter k.
    """
    n = len(v)
    
    # c_p counts how many T_i matrices are e_p e_p^T
    # c_0 counts how many T_i are the zero matrix
    c_p = np.zeros(n, dtype=int)
    c_0 = 0
    
    for T in T_matrices:
        non_zero_elements = np.nonzero(np.diag(T))[0]
        if len(non_zero_elements) == 0:
            c_0 += 1
        elif len(non_zero_elements) == 1:
            # Assuming valid input where non-zero is 1
            p = non_zero_elements[0]
            c_p[p] += 1
        else:
            # This case should not happen based on problem description
            raise ValueError("Matrix has more than one non-zero diagonal entry.")

    total_value = 0.0
    
    print("The final value is calculated by the formula:")
    print("Value = sum_{p=0}^{n-1} [ (c_p + c_0)^k - c_0^k ] * S_pp * v_p^2\n")
    print("Calculation breakdown:")
    
    full_equation = []
    
    for p in range(n):
        # Number of sequences where all matrices map to position p or are zero
        num_sequences_p_or_0 = (c_p[p] + c_0) ** k
        # Number of sequences where all matrices are zero
        num_sequences_0 = c_0 ** k
        
        # Number of sequences where the product results in e_p e_p^T
        count = num_sequences_p_or_0 - num_sequences_0
        
        # Value for each such sequence
        term_value = S[p, p] * (v[p] ** 2)
        
        # Total contribution from position p
        total_contribution = count * term_value
        total_value += total_contribution
        
        print(f"p={p}:")
        print(f"  c_{p} = {c_p[p]}, c_0 = {c_0}")
        print(f"  Sequences count = ({c_p[p]} + {c_0})^{k} - {c_0}^{k} = {count}")
        print(f"  Term value = S[{p},{p}] * v[{p}]^2 = {S[p,p]} * {v[p]**2} = {term_value}")
        print(f"  Contribution = {count} * {term_value} = {total_contribution}")
        
        # Storing parts of the final equation
        full_equation.append(f"{total_contribution}")

    print("\nFinal Equation:")
    # Printing the full summation
    final_eq_str = " + ".join(full_equation)
    print(f"Value = {final_eq_str} = {total_value}")

    print(f"\nFinal calculated VEST value: {total_value}")


# Example Usage:
n_dim = 3
k_param = 2

# v in Q^n
v_vec = np.array([2.0, -1.0, 3.0])

# T_i are diagonal Z_2 matrices with at most one non-zero entry.
# c_0=1 (T1), c_1=2 (T2, T4), c_2=1 (T3)
T_list = [
    np.diag(np.array([1.0, 0, 0])),  # Corresponds to p=0
    np.diag(np.array([0, 1.0, 0])),  # Corresponds to p=1
    np.diag(np.array([0, 0, 1.0])),  # Corresponds to p=2
    np.diag(np.array([0, 1.0, 0])),  # Corresponds to p=1 again
    np.diag(np.array([0, 0, 0]))    # The zero matrix
]

# S in Q^{n x n}
S_mat = np.array([
    [1.0, 0.5, 0],
    [0.5, 2.0, 0],
    [0,   0,   0.5]
])

solve_vest_b(v_vec, T_list, S_mat, k_param)