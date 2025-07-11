import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def main():
    """
    Computes the order of the Galois group for f(x) = x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # This is a depressed quartic of the form x^4 + px + q.
    p = 8
    q = 14

    print("Computing the order of the Galois group for the polynomial f(x) = x^4 + 8x + 14.")
    print("---")

    # Step 1: Check for irreducibility over Q
    print("Step 1: Check if the polynomial is irreducible over Q.")
    print("Using Eisenstein's criterion with the prime p = 2:")
    print("1. The prime 2 divides the coefficients 8 and 14.")
    print("2. The prime 2 does not divide the leading coefficient 1.")
    print("3. The prime's square, 4, does not divide the constant term 14.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("The Galois group is a transitive subgroup of S4: S4, A4, D4, V4, or C4.")
    print("---")

    # Step 2: Compute the discriminant
    print("Step 2: Compute the discriminant Δ.")
    print("For a polynomial x^4 + px + q, the discriminant is Δ = 256 * q^3 - 27 * p^4.")
    delta_val = 256 * (q**3) - 27 * (p**4)
    print(f"Δ = 256 * ({q}^3) - 27 * ({p}^4)")
    print(f"Δ = 256 * {q**3} - 27 * {p**4}")
    print(f"Δ = {256 * (q**3)} - {27 * (p**4)}")
    print(f"Δ = {delta_val}")
    print()

    # Check if the discriminant is a perfect square
    print(f"Checking if Δ = {delta_val} is a perfect square...")
    if is_perfect_square(delta_val):
        sqrt_delta = int(math.sqrt(delta_val))
        print(f"Yes, {delta_val} is a perfect square ({sqrt_delta}^2).")
        print("This means the Galois group is a subgroup of A4 (A4 or V4).")
    else:
        print(f"No, {delta_val} is not a perfect square (sqrt({delta_val}) is approx {math.sqrt(delta_val):.2f}).")
        print("This means the Galois group is not a subgroup of A4. Possible groups are S4 or D4.")
    print("---")

    # Step 3: Construct and analyze the resolvent cubic
    print("Step 3: Analyze the resolvent cubic g(y).")
    # For f(x) = x^4 + ax^3 + bx^2 + cx + d
    # Here, a=0, b=0, c=8, d=14.
    # The resolvent cubic is y^3 - by^2 + (ac-4d)y - (a^2d - 4bd + c^2).
    a, b, c, d = 0, 0, 8, 14
    resolvent_coeff_y2 = -b
    resolvent_coeff_y1 = a * c - 4 * d
    resolvent_coeff_y0 = -(a**2 * d - 4 * b * d + c**2)
    print(f"For f(x), a={a}, b={b}, c={c}, d={d}. The resolvent cubic is calculated as:")
    print(f"g(y) = y^3 - ({b})y^2 + (({a})*({c}) - 4*({d}))y - (({a}^2)*({d}) - 4*({b})*({d}) + {c}^2)")
    print(f"g(y) = y^3 - ({resolvent_coeff_y2})y^2 + ({resolvent_coeff_y1})y + ({resolvent_coeff_y0})")
    print("So, the resolvent cubic is g(y) = y^3 - 56y - 64.")
    print()

    # Check if the resolvent cubic is reducible over Q
    print("Checking if g(y) is reducible over Q by finding rational roots.")
    print("By the Rational Root Theorem, any rational root must be an integer divisor of -64.")
    # The constant term is -64. Divisors are +/-1, +/-2, +/-4, +/-8, etc. Let's test them.
    root_candidate = 8
    val_at_candidate = root_candidate**3 + resolvent_coeff_y1 * root_candidate + resolvent_coeff_y0
    
    print(f"Testing a divisor, y = {root_candidate}:")
    print(f"g({root_candidate}) = ({root_candidate})^3 - 56*({root_candidate}) - 64")
    print(f"g({root_candidate}) = {root_candidate**3} - {56*root_candidate} - 64")
    print(f"g({root_candidate}) = {val_at_candidate}")
    if val_at_candidate == 0:
        print("Since g(8) = 0, the resolvent cubic is reducible over Q.")
    else:
        # This part of the code won't be reached for this specific polynomial
        print("The resolvent cubic is irreducible over Q.")
    print("---")

    # Step 4: Final Conclusion
    print("Step 4: Determine the Galois group.")
    print("Summary of findings:")
    print("1. f(x) is irreducible over Q.")
    print(f"2. The discriminant Δ = {delta_val} is NOT a perfect square.")
    print("3. The resolvent cubic g(y) IS reducible over Q.")
    print("These conditions uniquely identify the Galois group of the polynomial as the Dihedral group D4.")
    print()
    print("The order of the Dihedral group D4 (the group of symmetries of a square) is 8.")

if __name__ == "__main__":
    main()
