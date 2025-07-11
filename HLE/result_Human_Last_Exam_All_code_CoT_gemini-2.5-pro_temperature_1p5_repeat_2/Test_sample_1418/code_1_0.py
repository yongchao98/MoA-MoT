import math

def solve_and_check():
    """
    This function solves the equation z^2 - 36z + 288 = 0
    and checks if the solutions for x, where z=x^2, are rational.
    """
    print("To check if beta = (2+sqrt(2))(3+sqrt(3)) is a square in Q(sqrt(2),sqrt(3)),")
    print("we check if its norm, 36 + 24*sqrt(2), is a square in Q(sqrt(2)).")
    print("This leads to the equation x^4 - 36*x^2 + 288 = 0 for a rational number x.")
    print("Let z = x^2, the equation becomes a quadratic in z:")
    
    a = 1
    b = -36
    c = 288

    print(f"z^2 + ({b})z + ({c}) = 0")
    
    # Discriminant of the quadratic equation
    discriminant = b**2 - 4*a*c
    print(f"The discriminant is ({b})^2 - 4*({a})*({c}) = {discriminant}")

    if discriminant < 0:
        print("The discriminant is negative, so there are no real solutions for z.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    print(f"The square root of the discriminant is {sqrt_discriminant:.2f}")
    
    z1 = (-b + sqrt_discriminant) / (2*a)
    z2 = (-b - sqrt_discriminant) / (2*a)
    print(f"The solutions for z = x^2 are z1 = {z1} and z2 = {z2}")

    # To find x, we take the square root of z.
    # For x to be rational, z must be a perfect square of a rational number.
    # We check if z1 or z2 are squares of rational numbers.
    def is_square_of_rational(n):
        if n < 0: return False
        # A rational number's square is rational. A rational number p/q has square p^2/q^2.
        # An integer is a perfect square if its square root is an integer.
        # Here z1=24, z2=12, are not squares of integers. And since they are integers,
        # their square roots are irrational.
        s = int(math.sqrt(n))
        return s*s == n

    if is_square_of_rational(z1) or is_square_of_rational(z2):
         print(f"At least one of the solutions for z, {z1} or {z2}, is a perfect square of a rational number.")
    else:
         print(f"Neither {z1} nor {z2} is the square of a rational number.")
         print("Therefore, there is no rational solution for x.")

    print("\nConclusion: This confirms that beta is not a square in Q(sqrt(2),sqrt(3)).")
    print("Thus, the degree of the extension L/Q is 8.")

solve_and_check()