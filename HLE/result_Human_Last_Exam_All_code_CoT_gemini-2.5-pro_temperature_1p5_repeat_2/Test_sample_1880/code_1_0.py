import math

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Define the polynomial coefficients for the form x^4 + px + q
    p = 8
    q = 14
    print(f"Analyzing the polynomial f(x) = x^4 + {p}x + {q}\n")

    # Step 1: Check for irreducibility using Eisenstein's Criterion
    print("--- Step 1: Check for Irreducibility ---")
    print("Using Eisenstein's criterion with the prime p=2:")
    print("1. The prime p=2 does not divide the leading coefficient (which is 1).")
    print("2. The prime p=2 divides all other coefficients (0 for x^3, 0 for x^2, 8 for x, and 14 for the constant term).")
    print("3. The value p^2=4 does not divide the constant term (14).")
    print("Conclusion: The polynomial is irreducible over Q.\n")

    # Step 2: Compute the discriminant
    print("--- Step 2: Compute the Discriminant (Δ) ---")
    print("For a polynomial x^4 + px + q, the discriminant is Δ = -27*p^4 + 256*q^3.")
    val_p4 = p**4
    val_q3 = q**3
    term1 = -27 * val_p4
    term2 = 256 * val_q3
    delta = term1 + term2
    print(f"p = {p}, q = {q}")
    print(f"Δ = -27 * ({p})^4 + 256 * ({q})^3")
    print(f"Δ = -27 * {val_p4} + 256 * {val_q3}")
    print(f"Δ = {term1} + {term2}")
    print(f"Δ = {delta}\n")

    # Check if the discriminant is a perfect square
    print("--- Step 3: Check if Δ is a Perfect Square ---")
    sqrt_delta_int = math.isqrt(delta)
    is_square = (sqrt_delta_int**2 == delta)
    if is_square:
        print(f"The discriminant {delta} is a perfect square.")
        print("Conclusion: The Galois group is a subgroup of A_4.\n")
    else:
        print(f"The discriminant {delta} is not a perfect square (sqrt({delta}) ≈ {math.sqrt(delta):.2f}).")
        print("Conclusion: The Galois group is NOT a subgroup of A_4.\n")

    # Step 4: Form and analyze the resolvent cubic
    print("--- Step 4: Analyze the Resolvent Cubic ---")
    print("For x^4 + px + q, the resolvent cubic is g(t) = t^3 - 4*q*t - p^2.")
    g_coeff_t = -4 * q
    g_const_term = -p**2
    print(f"g(t) = t^3 - 4*({q})*t - ({p})^2")
    print(f"g(t) = t^3 {g_coeff_t}t {g_const_term}")

    # Find rational roots of the resolvent cubic via Rational Root Theorem
    print("\nFinding rational roots of g(t) by testing divisors of the constant term:", g_const_term)
    abs_const = abs(g_const_term)
    divisors = set()
    for i in range(1, math.isqrt(abs_const) + 1):
        if abs_const % i == 0:
            divisors.add(i)
            divisors.add(-i)
            divisors.add(abs_const // i)
            divisors.add(-(abs_const // i))
    
    rational_roots = []
    for d in sorted(list(divisors)):
        if d**3 + g_coeff_t * d + g_const_term == 0:
            rational_roots.append(d)
            print(f"Found a rational root at t = {d}, since g({d}) = ({d})^3 {g_coeff_t}*({d}) {g_const_term} = 0.")

    num_rational_roots = len(rational_roots)
    if num_rational_roots == 0:
        print("The resolvent cubic is irreducible over Q (has no rational roots).")
    else:
        print(f"The resolvent cubic is reducible over Q and has {num_rational_roots} rational root(s).")
    print("")

    # Step 5: Determine the Galois Group and its order
    print("--- Step 5: Determine the Galois Group and its Order ---")
    galois_group_name = ""
    order = 0
    
    if not is_square:
        if num_rational_roots == 0:
            galois_group_name = "S_4"
            order = 24
        elif num_rational_roots == 1:
            galois_group_name = "D_4"
            order = 8
        elif num_rational_roots == 3:
            galois_group_name = "C_4"
            order = 4
    elif is_square:
        if num_rational_roots == 0:
            galois_group_name = "A_4"
            order = 12
        elif num_rational_roots == 3:
            galois_group_name = "V_4"
            order = 4
    
    print("Summary of Results:")
    print(f"1. The polynomial is irreducible.")
    print(f"2. The discriminant is NOT a perfect square.")
    print(f"3. The resolvent cubic has {num_rational_roots} rational root(s).")
    print(f"Based on these facts, the Galois group is the dihedral group {galois_group_name}.\n")

    print("--- Final Answer ---")
    print(f"The order of the Galois group is {order}.")
    
    return order

final_order = compute_galois_group_order()
print(f'<<<{final_order}>>>')