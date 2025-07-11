def N(a):
    """
    Calculates the numerator of the simple continued fraction [a_1, a_2, ...].
    The recurrence relation is p_n = a_n * p_{n-1} + p_{n-2}
    with initial conditions p_0 = 1, p_{-1} = 0.
    """
    p_nm2 = 0
    p_nm1 = 1
    for x in a:
        p_n = x * p_nm1 + p_nm2
        p_nm2 = p_nm1
        p_nm1 = p_n
    return p_nm1

def solve_and_verify(a):
    """
    Solves for c_k and verifies the identity for a given sequence a = [a_1, ..., a_k].
    """
    k = len(a)
    if k < 2:
        print("k must be >= 2.")
        return

    print(f"Verifying for k={k} and a = {a}")
    
    # Construct the sequences for the LHS and RHS
    a_LHS = a[1:-1] + [a[-1] + 1] + list(reversed(a))
    a_RHS = a + list(reversed(a))[1:]

    # Calculate the values from the equation
    val_LHS = N(a_LHS)
    val_RHS = N(a_RHS)
    
    # Calculate c_k from the original equation
    c_k_from_eq = val_LHS - val_RHS

    # Calculate c_k from the derived formula
    # c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]
    term1 = N(a)
    # For k=2, a[1:-1] is empty, N([])=1, which is correct for q_{k-1}
    term2 = N(a[1:-1]) if k > 2 else N([])
    c_k_from_formula = term1 * term2

    # Print the results
    print(f"LHS: N({a_LHS}) = {val_LHS}")
    print(f"RHS term: N({a_RHS}) = {val_RHS}")
    print(f"From the equation, c_{k} = {val_LHS} - {val_RHS} = {c_k_from_eq}")
    print("-" * 20)
    print(f"From the derived formula, c_{k} = N({a}) * N({a[1:-1] if k > 2 else '[]'})")
    print(f"c_{k} = {term1} * {term2} = {c_k_from_formula}")
    print("-" * 20)
    
    if c_k_from_eq == c_k_from_formula:
        print("The formula is verified successfully.")
        print("The final equation is:")
        print(f"{val_LHS} = {val_RHS} + {c_k_from_formula}")
    else:
        print("There is a discrepancy in the verification.")

# Example with k=3, a = [2, 3, 4]
a_example = [2, 3, 4]
solve_and_verify(a_example)

print("\n" + "="*40 + "\n")

# Example with k=2, a = [5, 6]
a_example_2 = [5, 6]
solve_and_verify(a_example_2)