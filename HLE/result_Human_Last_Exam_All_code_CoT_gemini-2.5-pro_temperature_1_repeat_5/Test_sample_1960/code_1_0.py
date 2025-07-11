def print_general_encoding_formulas():
    """
    This function prints the general formulas for f(w) and C(m, b) that encode
    the equipartitioning problem in multiplicative linear logic.
    
    Notation:
    _|_ : bottom constant
    --o : linear implication (multimap)
    @   : tensor product
    ^n  : formula tensored with itself n times (e.g., A^2 is A @ A)
    """

    m = "m"  # symbolic variable for the number of partitions
    b = "b"  # symbolic variable for the target sum
    w = "w"  # symbolic variable for a number in the set W
    
    # 1. Define the base formula A
    A = "(_|_ --o _|_)"
    
    # 2. Define the formula for f(w)
    # f(w) = A^w
    f_w_symbolic = f"({A})^({w})"
    
    # 3. Define the formula for C(m, b)
    # It is built from a verifier formula P_b = ((A^b) --o _|_) --o _|_
    A_pow_b_symbolic = f"({A})^({b})"
    P_b_symbolic = f"((({A_pow_b_symbolic}) --o _|_) --o _|_)"
    # C(m,b) = (P_b)^m
    C_m_b_symbolic = f"({P_b_symbolic})^({m})"

    print("The encoding of the equipartitioning problem EP(W, m, b) is as follows:")
    print("--------------------------------------------------------------------")
    
    print("\nLet the base formula A be:")
    print(f"A = {A}")
    
    print(f"\nThe function f(w) that maps a number w to a formula is:")
    print(f"f({w}) = A^{w} = {f_w_symbolic}")
    print("(For w=0, f(0) = 1, the multiplicative unit)")
    
    print(f"\nThe goal formula C(m, b) for m={m} partitions and a target sum of b={b} is:")
    print(f"C({m}, {b}) = (P_b)^{m}, where P_b is a verifier formula.")
    print("\nThe verifier P_b is defined as:")
    print(f"P_b = ((A^b) --o _|_) --o _|_ = {P_b_symbolic}")
    print("\nSo, the complete formula for C is:")
    print(f"C({m}, {b}) = {C_m_b_symbolic}")
    print("--------------------------------------------------------------------")
    print("\nThe problem EP(W, m, b) is true iff the sequent {f(w) | w in W} |-- C(m, b) is derivable.")


print_general_encoding_formulas()