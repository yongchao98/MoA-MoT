import math

def N(a):
    """
    Calculates the numerator of the continued fraction defined by the list a.
    """
    if not a:
        return 1
    if len(a) == 1:
        return a[0]
    
    p_prev2 = 1
    p_prev1 = a[0]
    
    for i in range(1, len(a)):
        p_curr = a[i] * p_prev1 + p_prev2
        p_prev2 = p_prev1
        p_prev1 = p_curr
        
    return p_prev1

def solve_for_ck_example(k, a):
    """
    Calculates all terms in the equation for a given k and list a,
    then solves for c_k.
    """
    if len(a) != k:
        print(f"Error: The length of list 'a' must be equal to k ({k}).")
        return

    # Construct the argument list for the LHS numerator
    # LHS: N[a_2, ..., a_{k}+1, a_k, ..., a_1]
    lhs_args = a[1:k]
    lhs_args[-1] += 1
    lhs_args.extend(reversed(a[0:k]))

    # Construct the argument list for the RHS numerator
    # RHS: N[a_1, ..., a_k, a_k, ..., a_2]
    rhs_args = a[0:k]
    rhs_args.extend(reversed(a[1:k]))
    
    # Calculate the numerators
    lhs_val = N(lhs_args)
    rhs_val = N(rhs_args)
    
    # Calculate c_k
    ck = lhs_val - rhs_val
    
    # Verify with the derived formula
    # ck_formula = N(a) * N(a[1:k-1])
    # print(f"Verification using derived formula: c_k = {ck_formula}")

    print(f"For k = {k} and a = {a}:")
    print(f"The equation is: N({lhs_args}) = N({rhs_args}) + c_{k}")
    print(f"Numerically: {lhs_val} = {rhs_val} + {ck}")
    
# Example calculation for k=3 and a=[1, 2, 3]
k = 3
a = [1, 2, 3]
solve_for_ck_example(k, a)
