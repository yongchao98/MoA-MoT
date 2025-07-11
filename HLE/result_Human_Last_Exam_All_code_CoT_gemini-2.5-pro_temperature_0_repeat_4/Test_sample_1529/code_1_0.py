import sympy

def solve_discriminant():
    """
    Calculates the discriminant of the polynomial f(x) = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    # Define the symbolic variable and the polynomial
    x = sympy.Symbol('x')
    f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

    print(f"The curve is defined by the equation: y^2 = {f}")
    print("The task is to find the discriminant of the polynomial f(x).")

    # The polynomial has a root at x = -1, so (x + 1) is a factor.
    # We can factor f(x) into g(x) * h(x).
    g = x + 1
    try:
        h, rem = sympy.div(f, g, x)
        if rem != 0:
            raise ValueError("Division had a non-zero remainder.")
    except ValueError as e:
        print(f"Error during factorization: {e}")
        # Fallback to direct calculation if factorization fails
        disc_f = sympy.discriminant(f, x)
        print(f"\nDirect calculation of the discriminant of f(x) is: {disc_f}")
        return

    print(f"\nThe polynomial f(x) can be factored as f(x) = g(x) * h(x), where:")
    print(f"g(x) = {g}")
    print(f"h(x) = {h}")

    # The discriminant of a product of polynomials is given by the formula:
    # disc(f) = disc(g) * disc(h) * Res(g, h)^2
    # We will now calculate each component.

    disc_g = sympy.discriminant(g, x)
    disc_h = sympy.discriminant(h, x)
    res_gh = sympy.resultant(g, h, x)

    # Calculate the final discriminant using the formula
    final_discriminant = disc_g * disc_h * res_gh**2

    print("\nWe use the formula: disc(f) = disc(g) * disc(h) * Resultant(g, h)^2")
    print(f"The discriminant of g(x) is: {disc_g}")
    print(f"The discriminant of h(x) is: {disc_h}")
    print(f"The resultant of g(x) and h(x) is: {res_gh}")

    # Output the final equation with all the numbers
    print(f"\nThe final equation for the discriminant is:")
    print(f"{final_discriminant} = {disc_g} * {disc_h} * ({res_gh})^2")
    
    print(f"\nThus, the minimal discriminant of the curve is {final_discriminant}.")

    # For verification, we can also compute the discriminant of f(x) directly.
    direct_disc_f = sympy.discriminant(f, x)
    print(f"\nVerification by direct calculation gives: {direct_disc_f}")

if __name__ == '__main__':
    solve_discriminant()
