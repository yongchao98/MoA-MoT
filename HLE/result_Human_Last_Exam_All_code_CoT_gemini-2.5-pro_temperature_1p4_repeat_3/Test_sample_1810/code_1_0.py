import sys

def numerator(arr):
    """
    Computes the numerator of a continued fraction given by the list arr.
    N[a_1, ..., a_k]
    """
    if not arr:
        return 1
    if len(arr) == 1:
        return arr[0]
    
    # Using recurrence p_k = a_k * p_{k-1} + p_{k-2}
    p_prev2 = 1
    p_prev1 = arr[0]
    for i in range(1, len(arr)):
        p_curr = arr[i] * p_prev1 + p_prev2
        p_prev2 = p_prev1
        p_prev1 = p_curr
    return p_prev1

def solve_ck(a):
    """
    Solves for c_k for a given sequence a = [a_1, ..., a_k].
    """
    k = len(a)
    if k < 2:
        print("k must be >= 2")
        return None

    # Calculate p_k, p_{k-1}, q_k, q_{k-1}
    p_k = numerator(a)
    p_k_minus_1 = numerator(a[:-1])
    
    q_k = numerator(a[1:])
    if k > 2:
        q_k_minus_1 = numerator(a[1:-1])
    elif k == 2:
        q_k_minus_1 = 1 # N[] for the sequence a[1:-1] which is empty
    
    ak = a[-1]

    # Expression for LHS: N[a_2,..., a_{k}+1, a_k,...,a_1]
    lhs_seq = a[1:] + [ak + 1] + list(reversed(a))
    lhs = numerator(lhs_seq)

    # Expression for RHS: N[a_1,...,a_{k}, a_k,...,a_2]
    rhs_seq = a + list(reversed(a[:-1]))
    rhs = numerator(rhs_seq)
    
    c_k = lhs - rhs
    
    # Print the values to verify the formula
    print(f"For the sequence a = {a}:")
    print(f"k = {k}")
    print(f"a_k = {ak}")
    print(f"p_k = N{a} = {p_k}")
    print(f"p_(k-1) = N{a[:-1]} = {p_k_minus_1}")
    print(f"q_k = N{a[1:]} = {q_k}")
    print(f"q_(k-1) = N{a[1:-1] if k > 2 else '[]'} = {q_k_minus_1}")
    
    # Calculate c_k using the derived formula
    # c_k = a_k*p_k*q_k + p_k*q_{k-1} + p_{k-1}*q_k - p_{k-1}*q_{k-1}
    formula_ck = ak * p_k * q_k + p_k * q_k_minus_1 + p_k_minus_1 * q_k - p_k_minus_1 * q_k_minus_1

    print(f"\nEquation:")
    print(f"N{[round(n,2) for n in lhs_seq]} = {lhs}")
    print(f"N{[round(n,2) for n in rhs_seq]} = {rhs}")
    print(f"c_k = LHS - RHS = {lhs} - {rhs} = {c_k}")
    print(f"Value from derived formula: {formula_ck}")

    if c_k == formula_ck:
        print("\nFormula is verified.")
        print("The expression for c_k is:")
        print("c_k = a_k * p_k * q_k + p_k * q_{k-1} + p_{k-1} * q_k - p_{k-1} * q_{k-1}")
        print("where:")
        print("p_k = N[a_1, ..., a_k]")
        print("p_{k-1} = N[a_1, ..., a_{k-1}]")
        print("q_k = N[a_2, ..., a_k]")
        print("q_{k-1} = N[a_2, ..., a_{k-1}]")

# Example usage with a sample list for a_1, ..., a_k
# You can change this list to test other sequences
sample_a = [1, 2, 3] 
solve_ck(sample_a)
