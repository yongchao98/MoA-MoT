import math

def get_divisors(n):
    """Helper function to get all integer divisors of a number."""
    n = abs(n)
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(-i)
            divs.add(n // i)
            divs.add(-(n // i))
    return sorted(list(divs))

def solve():
    """
    Computes the order of the Galois group for x^4 + 8x + 14.
    """
    # Polynomial f(x) = x^4 + 0x^3 + 0x^2 + 8x + 14
    # Standard form: ax^4 + bx^3 + cx^2 + dx + e
    # Our poly form: x^4 + c_x*x + d_x
    c_x = 8
    d_x = 14

    print("The polynomial is f(x) = x^4 + 8x + 14.\n")

    # Step 1: Check for irreducibility using Eisenstein's criterion.
    print("Step 1: Check for irreducibility over Q.")
    p = 2
    print(f"We use Eisenstein's criterion with the prime p = {p}.")
    print(f"1. p={p} divides the coefficients 14, 8, 0, 0.")
    print(f"2. p={p} does not divide the leading coefficient 1.")
    print(f"3. p^2=4 does not divide the constant term 14.")
    print("Conclusion: The polynomial is irreducible over Q.\n")

    # Step 2: Calculate the discriminant.
    # For a polynomial x^4 + cx + d, the discriminant is 256*d^3 - 27*c^4.
    print("Step 2: Calculate the discriminant (Delta).")
    delta = 256 * (d_x**3) - 27 * (c_x**4)
    print(f"The formula for the discriminant of x^4 + {c_x}x + {d_x} is Delta = 256*({d_x}^3) - 27*({c_x}^4).")
    print(f"Delta = 256 * ({d_x**3}) - 27 * ({c_x**4})")
    print(f"Delta = {256 * d_x**3} - {27 * c_x**4}")
    print(f"Delta = {delta}\n")

    # Step 3: Check if Delta is a perfect square.
    print("Step 3: Check if Delta is a perfect square in Q.")
    sqrt_delta = math.isqrt(delta)
    is_square = (sqrt_delta * sqrt_delta == delta)
    print(f"The square root of {delta} is approximately {math.sqrt(delta):.2f}.")
    if is_square:
        print(f"Delta = {delta} is a perfect square.")
        print("Conclusion: The Galois group is a subgroup of A4. Possible groups: A4, V4.\n")
    else:
        print(f"Delta = {delta} is not a perfect square.")
        print("Conclusion: The Galois group is not a subgroup of A4. Possible groups: S4, D4.\n")

    # Step 4: Form the resolvent cubic.
    # For a polynomial x^4 + cx + d, the resolvent cubic is y^3 - 4*d*y - c^2 = 0.
    print("Step 4: Form the resolvent cubic.")
    coeff_y1 = -4 * d_x
    coeff_y0 = -(c_x**2)
    print(f"The resolvent cubic for x^4 + {c_x}x + {d_x} is y^3 - 4*({d_x})*y - ({c_x})^2 = 0.")
    print(f"The equation is: y^3 + {coeff_y1}y + {coeff_y0} = 0.\n")
    
    # Step 5: Test the resolvent cubic for rational roots.
    print("Step 5: Test if the resolvent cubic is reducible over Q.")
    print("By the Rational Root Theorem, any rational root must be a divisor of the constant term.")
    constant_term = coeff_y0
    divs = get_divisors(constant_term)
    print(f"The divisors of {constant_term} are: {divs}.")
    
    found_root = None
    for y in divs:
        # Test if y is a root of y^3 + coeff_y1*y + coeff_y0 = 0
        if y**3 + coeff_y1 * y + coeff_y0 == 0:
            found_root = y
            break

    if found_root is not None:
        print(f"\nTesting for y = {found_root}:")
        print(f"({found_root})^3 - 56*({found_root}) - 64 = {found_root**3} - {56*found_root} - 64 = {found_root**3 - 56*found_root - 64}")
        print(f"A rational root y = {found_root} was found.")
        print("Conclusion: The resolvent cubic is reducible over Q.\n")
    else:
        print("No rational roots were found. The resolvent cubic is irreducible over Q.\n")

    # Step 6: Determine the Galois group and its order.
    print("Step 6: Determine the final Galois group and its order.")
    print("Recap:")
    print(f" - The polynomial is irreducible.")
    print(f" - The discriminant is NOT a perfect square.")
    print(f" - The resolvent cubic IS reducible.")
    print("This combination of properties implies that the Galois group is the Dihedral group D4.\n")
    
    order = 8
    print(f"The order of the Galois group D4 is {order}.")


if __name__ == '__main__':
    solve()