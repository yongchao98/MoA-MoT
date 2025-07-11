import sympy

def solve_curve_thickness():
    """
    Calculates the 2-adic valuation of the discriminant of the polynomial
    defining the curve z^2 = 2*x^5 + 2*x^3 + 1.
    This valuation corresponds to the thickness of the double point in the stable reduction.
    """
    # Define the symbolic variable x
    x = sympy.Symbol('x')

    # Define the polynomial f(x)
    f = 2 * x**5 + 2 * x**3 + 1

    # Calculate the discriminant of the polynomial
    poly_discriminant = sympy.discriminant(f, x)

    # The prime for the valuation is 2
    p = 2
    
    # Calculate the p-adic valuation of the discriminant
    # The valuation is the exponent of p in the prime factorization of the discriminant.
    # We can find it by repeatedly dividing the discriminant by p until it's no longer divisible.
    val = 0
    num = int(poly_discriminant)
    
    # We need to handle the case where the discriminant is zero or negative.
    if num == 0:
        # A discriminant of 0 means the polynomial has repeated roots,
        # but the valuation would be infinite. We print the calculation.
        print(f"The polynomial is f(x) = {f}")
        print(f"The discriminant is {num}")
        print("The discriminant is 0, so the valuation is infinite.")
        return

    temp_num = abs(num)
    while temp_num > 0 and temp_num % p == 0:
        val += 1
        temp_num //= p
        
    print(f"The curve is defined by the equation: z^2 = 2*x^5 + 2*x^3 + 1")
    print(f"The polynomial is f(x) = {f}")
    print(f"The discriminant of the polynomial f(x) is: {poly_discriminant}")
    print(f"The prime factorization of the discriminant is: 2^{val} * {num // (p**val)}")
    print(f"The 2-adic valuation of the discriminant is the exponent of 2, which is: {val}")
    print(f"The thickness of the double point of the stable reduction above 2 is equal to this valuation.")
    print(f"Therefore, the thickness is {val}.")

solve_curve_thickness()