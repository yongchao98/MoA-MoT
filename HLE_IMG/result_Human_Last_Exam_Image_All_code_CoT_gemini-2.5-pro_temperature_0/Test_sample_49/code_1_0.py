import sympy as sp

def solve_cutoff_frequency():
    """
    This function symbolically derives the cutoff frequency for the given ladder network.
    """
    # Define symbols for resistance, capacitance, and the ladder's equivalent resistance
    r, C, R_ladder = sp.symbols('r C R_ladder', positive=True, real=True)

    # --- Step 1: Find the equivalent resistance of the infinite ladder (R_ladder) ---
    # The ladder's resistance is described by the recursive formula:
    # R_ladder = r (series) + r (series) + (r (shunt) || R_ladder (next section))
    # R_ladder = 2*r + (r * R_ladder) / (r + R_ladder)
    # This leads to a quadratic equation for R_ladder:
    # R_ladder^2 - 2*r*R_ladder - 2*r^2 = 0
    equation = sp.Eq(R_ladder**2 - 2*r*R_ladder - 2*r**2, 0)
    solutions = sp.solve(equation, R_ladder)

    # Resistance must be positive, so we take the positive root.
    # The solutions are r*(1 - sqrt(3)) and r*(1 + sqrt(3)).
    R_ladder_sol = solutions[1]

    # --- Step 2: Find the Thevenin resistance (R_th) at node a0 ---
    # To find R_th, the input source is grounded (c0 is at ground potential).
    # R_th is the resistance from a0 to ground. It consists of one resistor 'r'
    # (from a0 to c0) in parallel with the rest of the network seen from a0.
    #
    # The resistance of the other path from a0 is: r (from a0 to d1) in series
    # with the resistance from d1 to ground (R_d1_gnd).
    #
    # Resistance from d1 to ground (R_d1_gnd):
    # This path goes from d1 to c1, then from c1 to ground.
    # The resistance between d1 and c1 is the first shunt resistor 'r' in parallel with R_ladder.
    R_d1c1 = (r * R_ladder_sol) / (r + R_ladder_sol)
    # The resistance from c1 to ground is 'r'.
    R_d1_gnd = R_d1c1 + r

    # The total resistance of the path starting from a0 through d1 is:
    R_path_d1 = r + R_d1_gnd

    # R_th is 'r' (from a0 to ground) in parallel with R_path_d1
    R_th = (r * R_path_d1) / (r + R_path_d1)
    R_th_simplified = sp.simplify(R_th) # This evaluates to r*(sqrt(3) - 1)

    # --- Step 3: Determine the cutoff frequency (f_c) ---
    # The cutoff frequency for an RC low-pass filter is f_c = 1 / (2 * pi * R_th * C)
    
    print("The analysis of the circuit yields the following results:")
    print(f"1. Equivalent resistance of the infinite ladder: R_ladder = {R_ladder_sol}")
    print(f"2. Thevenin equivalent resistance at node a0: R_th = {R_th_simplified}")
    
    print("\n" + "="*55)
    print("The final equation for the cutoff frequency f_c is:")
    # We construct the final expression string for clarity
    final_expression_str = f"f_c = 1 / (2 * pi * ({R_th_simplified}) * C)"
    print(final_expression_str)
    print("="*55 + "\n")
    
    print("The numbers and symbols in this final equation are:")
    print("- 1: The numerator of the fraction.")
    print("- 2: A constant factor in the denominator.")
    print("- pi: The mathematical constant pi (approx. 3.14159).")
    print("- sqrt(3): The square root of 3 (a constant, approx. 1.732).")
    print("- 1: A constant being subtracted from sqrt(3).")
    print("- r: The resistance value of the individual resistors.")
    print("- C: The capacitance value of the capacitor.")

if __name__ == '__main__':
    solve_cutoff_frequency()