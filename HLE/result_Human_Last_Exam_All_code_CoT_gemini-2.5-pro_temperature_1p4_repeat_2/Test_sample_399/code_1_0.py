import sympy

def solve_ivp():
    """
    Solves the initial value problem x'(t) = (t-1)x^2(t) with x(0) = -8, and finds x(1).
    """
    # Define initial conditions and target time
    t0 = 0
    x0 = -8
    t_target = 1

    print("The initial value problem is x'(t) = (t-1) * x^2(t), with x(0) = -8.")
    print("This is a separable differential equation.\n")

    print("Step 1: Separate variables and integrate.")
    print("∫(1/x^2) dx = ∫(t-1) dt")
    print("This gives the general solution: -1/x = t^2/2 - t + C\n")

    print("Step 2: Use the initial condition x(0) = -8 to find the constant C.")
    # Calculate C
    # -1/x0 = t0^2/2 - t0 + C
    # C = -1/x0 - (t0^2/2 - t0)
    C = -1/x0 - (t0**2 / 2 - t0)
    print(f"Substitute t = {t0} and x = {x0}:")
    print(f"-1/({x0}) = {t0}^2/2 - {t0} + C")
    print(f"1/8 = 0 + C")
    print(f"So, C = {sympy.Rational(C)}\n")

    print("Step 3: Write the particular solution.")
    print(f"-1/x(t) = t^2/2 - t + {sympy.Rational(C)}")
    print("Solving for x(t):")
    print("x(t) = 1 / (-(t^2/2 - t + 1/8))")
    print("x(t) = 1 / (-t^2/2 + t - 1/8)")
    print("Multiplying numerator and denominator by 8 for a cleaner look:")
    print("x(t) = 8 / (-4*t^2 + 8*t - 1)\n")


    print(f"Step 4: Find x({t_target}).")
    # Calculate x(t_target)
    # Denominator = -t_target**2/2 + t_target - C
    denominator = -(t_target**2 / 2 - t_target + C)
    x_target = 1 / denominator

    # Using the cleaned-up formula for printing
    num_clean = 8
    den_clean_term1 = -4 * t_target**2
    den_clean_term2 = 8 * t_target
    den_clean_term3 = -1
    den_clean = den_clean_term1 + den_clean_term2 + den_clean_term3

    print(f"Substitute t = {t_target} into the particular solution:")
    print(f"x({t_target}) = 8 / (-4*({t_target})^2 + 8*({t_target}) - 1)")
    print("Let's calculate the denominator:")
    print(f"Denominator = {den_clean_term1} + {den_clean_term2} + {den_clean_term3}")
    print(f"Denominator = {den_clean}")
    print("\nSo the final equation is:")
    print(f"x({t_target}) = {num_clean} / {den_clean}")

    # Use sympy.Rational for exact fraction representation
    result_fraction = sympy.Rational(num_clean, den_clean)
    result_float = float(x_target)
    print(f"\nThe value of x(1) is {result_fraction} or approximately {result_float:.4f}.")

solve_ivp()