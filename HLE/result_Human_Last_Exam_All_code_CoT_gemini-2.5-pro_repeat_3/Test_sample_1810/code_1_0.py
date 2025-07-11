import sys

def get_convergents(coeffs):
    """
    Computes the numerators (p) and denominators (q) for a continued fraction.
    p_k = a_k * p_{k-1} + p_{k-2}
    q_k = a_k * q_{k-1} + q_{k-2}
    Initial conditions: p_0 = 1, p_{-1} = 0 and q_0 = 0, q_{-1} = 1 are used for recurrence.
    The convergents for [a_1, ..., a_k] are p_1/q_1, ..., p_k/q_k.
    The function returns two lists: [p_1, ..., p_k] and [q_1, ..., q_k].
    """
    p_minus_1, p_0 = 0, 1
    q_minus_1, q_0 = 1, 0
    
    p_list = []
    q_list = []
    
    for a in coeffs:
        p_k = a * p_0 + p_minus_1
        q_k = a * q_0 + q_minus_1
        
        p_list.append(p_k)
        q_list.append(q_k)
        
        p_minus_1, p_0 = p_0, p_k
        q_minus_1, q_0 = q_0, q_k
        
    return p_list, q_list

def get_numerator(coeffs):
    """
    Computes the numerator of a continued fraction defined by a list of coeffs.
    """
    if not coeffs:
        return 1
    
    p_list, _ = get_convergents(coeffs)
    return p_list[-1]

def solve_c_k(a):
    """
    Solves for c_k given a list of coefficients a = [a_1, ..., a_k].
    """
    k = len(a)
    if k < 2:
        print("Error: k must be >= 2.")
        return

    # Calculate convergents for [a_1, ..., a_k]
    p_list, q_list = get_convergents(a)
    
    # p_k is the last element of p_list
    p_k = p_list[-1]
    
    # q_{k-1} is the second to last element of q_list
    q_k_minus_1 = q_list[-2]
    
    # The solution for c_k is p_k * q_{k-1}
    c_k = p_k * q_k_minus_1

    # For verification, let's compute the LHS and RHS of the equation.
    
    # LHS sequence: [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    lhs_seq = []
    if k > 2:
      lhs_seq.extend(a[1:-1]) # a_2 to a_{k-1}
    lhs_seq.append(a[-1] + 1)   # a_k + 1
    lhs_seq.extend(reversed(a)) # a_k to a_1
    
    lhs_val = get_numerator(lhs_seq)

    # RHS sequence: [a_1, ..., a_k, a_k, ..., a_2]
    rhs_seq = []
    rhs_seq.extend(a)         # a_1 to a_k
    if k > 1:
        rhs_seq.extend(reversed(a[1:])) # a_k to a_2

    rhs_val = get_numerator(rhs_seq)

    print(f"For k = {k} and a = {a}:")
    print(f"The solution is c_{k} = p_k * q_{k-1}")
    print(f"p_k = N{a} = {p_k}")
    print(f"q_{k-1} = D{a[:-1]} = {q_k_minus_1}")
    print(f"c_{k} = {p_k} * {q_k_minus_1} = {c_k}")
    print("\nVerification:")
    print(f"LHS = N{lhs_seq} = {lhs_val}")
    print(f"RHS = N{rhs_seq} + c_{k} = {rhs_val} + {c_k} = {rhs_val + c_k}")
    if lhs_val == rhs_val + c_k:
        print("The equation holds true.")
        print(f"The equation is: {lhs_val} = {rhs_val} + {c_k}")
    else:
        print("There is a discrepancy in the formula verification.")


# --- Main execution ---
# You can change the list 'a' to test with different coefficients and k.
# The length of 'a' is k.
# Example for k=3, a = [a1, a2, a3] = [1, 2, 3]
try:
    # Read coefficients from command line arguments, e.g., python script.py 1 2 3
    a_coeffs = [int(arg) for arg in sys.argv[1:]]
    if not a_coeffs:
        a_coeffs = [1, 2, 3] # Default values
except (ValueError, IndexError):
    print("Usage: python your_script_name.py a1 a2 a3 ...")
    print("Using default values a = [1, 2, 3] for demonstration.")
    a_coeffs = [1, 2, 3]

solve_c_k(a_coeffs)
<<<N[a_1, ..., a_k] * D[a_1, ..., a_{k-1}]>>>