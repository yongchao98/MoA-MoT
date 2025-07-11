import math

def solve_for_l():
    """
    This function performs all the steps to calculate the value of l.
    """
    print("Step 1 & 2: Calculate u1 and the function terms.")
    
    # Given u2, calculate u1 from the center of gravity of the welded sheet.
    # The welded sheet consists of a rectangle B [0, 2a] x [0, a] and a rectangle C [0, 2a] x [a, 5a].
    # ks = (M_B*k_B + M_C*k_C) / (M_B + M_C)
    # 2a = (u1*2a^2 * a/2 + u2*8a^2 * 3a) / (u1*2a^2 + u2*8a^2)
    # This simplifies to 3*u1 = 8*u2
    u2 = 3
    u1 = (8 / 3) * u2
    print(f"Given u2 = {u2}, we find u1 = {u1:.1f}")

    # Define the integrand g(t), f'(x), and f''(x)
    def g(t):
        return (2 * t**3 + t) / (1 + t**4)

    def f_prime(x):
        return (2 * x**3 + x) / (1 + x**4)

    def f_double_prime(x):
        numerator = (6 * x**2 + 1) * (1 + x**4) - (2 * x**3 + x) * (4 * x**3)
        denominator = (1 + x**4)**2
        return numerator / denominator

    # Calculate f(5) using Simpson's rule with n=10
    n_simpson = 10
    a_bound, b_bound = 0, 5
    h = (b_bound - a_bound) / n_simpson
    integral = g(a_bound) + g(b_bound)
    for i in range(1, n_simpson, 2):
        integral += 4 * g(a_bound + i * h)
    for i in range(2, n_simpson, 2):
        integral += 2 * g(a_bound + i * h)
    integral *= h / 3
    f5_val = integral
    
    # Calculate f'(5) and f''(5)
    f_prime5_val = f_prime(5)
    f_double_prime5_val = f_double_prime(5)
    
    # Round the results to one decimal place
    f5_rounded = round(f5_val, 1)
    f_prime5_rounded = round(f_prime5_val, 1)
    f_double_prime5_rounded = round(f_double_prime5_val, 1)
    
    print(f"f(5) via Simpson's rule is ≈ {f5_val:.4f}, which rounds to {f5_rounded}")
    print(f"f'(5) is ≈ {f_prime5_val:.4f}, which rounds to {f_prime5_rounded}")
    print(f"f''(5) is ≈ {f_double_prime5_val:.4f}, which rounds to {f_double_prime5_rounded}")
    print("-" * 30)

    print("Step 3: Calculate the constant 'a'.")
    term = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded
    a_val = (u1 / 27) * (term**3)
    print(f"The expression (f(5) - 2f'(5) + 2f''(5)) with rounded values is: {f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}) = {term}")
    print(f"a = ({u1}/27) * ({term})^3 = {a_val}")
    print("-" * 30)

    print("Step 4 & 5: Calculate 'l'.")
    # For sheet A, the condition y_s = 4a simplifies to the relation:
    # 96 * a^3 = 2 * a * l^2
    # This further simplifies to l^2 = 48 * a^2, for a > 0.
    
    print("The final equation for l is derived from the center of gravity calculation for sheet A:")
    print("Equation: 96 * a^3 = 2 * a * l^2")
    
    print("\nSubstituting the calculated value of a:")
    lhs = 96 * a_val**3
    rhs_coeff = 2 * a_val
    print(f"96 * ({a_val:.1f})^3 = 2 * ({a_val:.1f}) * l^2")
    print(f"{lhs:.1f} = {rhs_coeff:.1f} * l^2")
    
    l_squared = lhs / rhs_coeff
    l_star = math.sqrt(l_squared)
    
    print(f"l^2 = {lhs:.1f} / {rhs_coeff:.1f} = {l_squared:.1f}")
    print(f"l = sqrt({l_squared:.1f}) = {l_star:.4f}")

    print("\nFinal Answer:")
    print(f"<<<{l_star}>>>")

# Run the solver
solve_for_l()