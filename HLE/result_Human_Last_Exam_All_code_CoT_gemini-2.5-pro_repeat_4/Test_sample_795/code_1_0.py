import math

def print_magnetization_formula():
    """
    Prints the derived analytical expression for the initial magnetization curve
    of a thin superconducting strip in a perpendicular magnetic field.
    """

    # --- Introduction to the result ---
    print("Based on the critical-state model for a thin superconducting strip (b << a) in a perpendicular magnetic field,")
    print("the analytical expression for the initial magnetization curve M(H) is derived.")
    print("\nThe final formula is:")

    # --- The Final Equation ---
    # The equation is: M(H) = - (Jc * a / 4) * tanh^2( (pi * H) / (2 * Jc * b) )
    # We will print it piece by piece to satisfy the "output each number" instruction.
    
    # Define the components of the equation symbolically
    term_Jc = "Jc"
    term_a = "a"
    term_b = "b"
    term_H = "H"
    term_pi = "π"
    term_tanh_sq = "[tanh(...)]^2"
    
    # Print the equation in a structured way
    print(f"\nM(H) = - ( ({term_Jc} * {term_a}) / 4 ) * {term_tanh_sq}")
    print(f"         where the argument of tanh is: ( ({term_pi} * {term_H}) / (2 * {term_Jc} * {term_b}) )")

    # --- Explanation of Terms ---
    print("\n--- Equation Components Explained ---")
    print("This equation describes the magnetization `M` as a function of the applied field `H`.")
    print("Here is a breakdown of each symbol and number in the final equation:")
    
    print("\nVariables and Constants:")
    print(f"  M(H): The Magnetization of the superconductor.")
    print(f"  {term_H}:   The strength of the uniform applied magnetic field.")
    print(f"  {term_Jc}:   The constant critical-current density of the superconductor.")
    print(f"  {term_a}:   The half-width of the bar's cross-section (dimension in x-direction).")
    print(f"  {term_b}:   The half-thickness of the bar's cross-section (dimension in y-direction).")
    print(f"  {term_pi}:   The mathematical constant pi (approximately {math.pi:.5f}).")

    print("\nNumerical Coefficients:")
    print(f"  4:    This number in the denominator `(Jc * a / 4)` comes from the geometric calculation of the magnetic moment.")
    print(f"  2:    This number in the denominator `(2 * Jc * b)` is part of the characteristic field that governs the penetration depth.")

    print("\nFunctions:")
    print(f"  tanh: The hyperbolic tangent function, which describes the saturation behavior.")
    print(f"  ^2:   The result of the tanh function is squared.")

if __name__ == "__main__":
    print_magnetization_formula()
    # The final expression in a single line for clarity.
    final_expression = "M(H) = -(Jc*a/4) * tanh^2((π*H)/(2*Jc*b))"
    print(f"\n<<<M(H) = -(Jc*a/4)*tanh^2((pi*H)/(2*Jc*b))>>>")