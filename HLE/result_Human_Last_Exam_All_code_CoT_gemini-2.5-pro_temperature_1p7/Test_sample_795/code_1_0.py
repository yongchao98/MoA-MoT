import math

def display_magnetization_formula():
    """
    Prints the derived analytical expression for the initial magnetization
    curve of a thin superconducting strip.
    """

    # --- Introduction and Symbol Definitions ---
    print("This script presents the analytical expression for the initial magnetization curve M(H) of a thin (b << a) superconducting strip in a perpendicular magnetic field H, based on the critical-state model.")
    print("-" * 70)
    print("The symbols used in the equation are defined as follows:")
    definitions = {
        "M(H)": "The magnetization of the superconductor.",
        "H": "The uniform applied magnetic field.",
        "Jc": "The constant critical-current density.",
        "a": "The half-width of the strip (along the x-axis).",
        "b": "The half-thickness of the strip (along the y-axis).",
        "tanh": "The hyperbolic tangent function.",
        "π": f"The mathematical constant pi (approx. {math.pi})."
    }
    for symbol, definition in definitions.items():
        print(f"  {symbol:<6} -> {definition}")
    print("-" * 70)

    # --- Final Equation ---
    print("\nThe analytical expression for the initial magnetization curve is:")
    final_equation = "M(H) = - (Jc * a / 2) * tanh^2( (π * H) / (2 * b * Jc) )"
    print(f"\n    {final_equation}\n")

    # --- Breaking down the components of the equation as requested ---
    print("The equation consists of the following components:")
    print("\n  1. The Saturation Magnetization Magnitude (M_sat):")
    print("     M_sat = Jc * a / 2")
    print("     This is the maximum possible magnetization when the current has penetrated the entire strip.")
    
    print("\n  2. The Argument of the Hyperbolic Tangent Function (x):")
    print("     x = (π * H) / (2 * b * Jc)")
    print("     This dimensionless term relates the applied field H to the characteristic field of the superconductor.")
    
    print("\nSo the equation can be seen as M(H) = -M_sat * tanh^2(x).")

if __name__ == "__main__":
    display_magnetization_formula()
