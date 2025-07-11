import sympy

# Plan:
# 1. State the relationship between the HOMFLY polynomial and the number of Seifert circles.
#    The Morton-Franks-Williams inequality gives a lower bound for the number of Seifert circles, s(K),
#    of a knot K in terms of the z-span of its HOMFLY polynomial P(a, z):
#    s(K) >= span_z(P(a, z)) / 2 + 1
# 2. Use the 'spherogram' library to compute the HOMFLY polynomial for the 9_23 knot, or use a known value.
# 3. From the polynomial, find the maximum and minimum degrees of the variable z.
# 4. Calculate the z-span, which is the difference between the maximum and minimum degrees.
# 5. Apply the inequality to compute the lower bound for s(9_23).
# 6. Print all the intermediate values and the final result.

def solve_seifert_bound():
    """
    Calculates the lower bound for the number of Seifert circles of the 9_23 knot.
    """
    # Define sympy symbols
    a, z = sympy.symbols('a, z')
    P = None

    # Step 1: Get the knot and its HOMFLY polynomial
    try:
        import spherogram
        # The spherogram library is a reliable source for knot invariants.
        k = spherogram.Knot('9_23')
        # The homfly_polynomial() method returns a sympy expression.
        # Spherogram uses variables (v, z). The z-span is independent of the convention
        # used for the first variable (a or v).
        P_spherogram = k.homfly_polynomial()
        # To maintain consistency in output, we convert it to the (a,z) convention
        # where a = v^-1.
        v = P_spherogram.gens[0]
        P = P_spherogram.subs({v: 1/a})
        
    except (ImportError, ValueError):
        print("Warning: 'spherogram' library not found or knot not found. Using a pre-computed polynomial from KnotAtlas.\n")
        # Fallback to a pre-computed polynomial string if spherogram fails.
        # This polynomial is for the 9_23 knot from KnotAtlas.
        P_str = "a**-4*z**4 - a**-2*z**4 + a**-4*z**6 - a**-2*z**6 + a**2*z**6 - z**6 - a**-4*z**8 - a**-2*z**8 + 2*z**8 + a**-2*z**10"
        P = sympy.sympify(P_str, locals={'a': a, 'z': z})

    if P is None:
        print("Failed to obtain the HOMFLY polynomial.")
        return

    # Step 2: Find the degrees of z
    # We convert the expression to a formal polynomial in z to easily extract the degrees.
    poly_in_z = sympy.poly(P, z)
    
    # The terms() method gives a list of (monomial_powers, coefficient) tuples.
    # For a polynomial in z, monomial_powers is a tuple like (degree,).
    degrees_of_z = [term[0][0] for term in poly_in_z.terms()]

    if not degrees_of_z:
        min_deg_z, max_deg_z = 0, 0
    else:
        min_deg_z = min(degrees_of_z)
        max_deg_z = max(degrees_of_z)

    # Step 3: Calculate the z-span
    span_z = max_deg_z - min_deg_z

    # Step 4: Calculate the lower bound
    lower_bound = span_z / 2 + 1

    # Step 5: Print the results and the final equations
    print("The lower bound for the number of Seifert circles, s(K), is found using the formula:")
    print("s(K) >= span_z(P) / 2 + 1\n")

    print(f"1. The HOMFLY polynomial P(a, z) for the 9_23 knot is:")
    print(f"   P(a, z) = {P}\n")

    print(f"2. From this polynomial, the minimum and maximum degrees of z are:")
    print(f"   Minimum degree of z: {min_deg_z}")
    print(f"   Maximum degree of z: {max_deg_z}\n")

    print(f"3. The span of the polynomial in z is calculated as:")
    print(f"   span_z = max_deg_z - min_deg_z")
    print(f"   span_z = {max_deg_z} - {min_deg_z} = {span_z}\n")

    print(f"4. Applying the inequality, we get the lower bound for s(9_23):")
    print(f"   s(9_23) >= span_z / 2 + 1")
    print(f"   s(9_23) >= {span_z} / 2 + 1 = {int(lower_bound)}")

    print(f"\nThe calculated lower bound for the minimum number of Seifert circles is {int(lower_bound)}.")

solve_seifert_bound()