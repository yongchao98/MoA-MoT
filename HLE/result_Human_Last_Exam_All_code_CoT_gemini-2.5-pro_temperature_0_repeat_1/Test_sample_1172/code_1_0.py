import math

def calculate_inductance_change():
    """
    Calculates the change in mutual inductance between two wire circuits
    when ideal magnetic concentrators are added.
    """

    # --- User-defined parameters ---
    # These values are chosen as an example satisfying d >> h.
    # Separation of wires in each circuit (meters)
    h = 0.01
    # Distance between the centers of the two circuits (meters)
    d = 0.1
    # Inner radius of the concentrator shells (meters)
    R1 = 0.02
    # Outer radius of the concentrator shells (meters)
    R2 = 0.04

    # --- Constants ---
    # Permeability of free space (Henry per meter)
    mu_0 = 4 * math.pi * 1e-7

    # --- Derivation and Final Expression ---
    print("This script calculates the change in mutual inductance (ΔM) per unit length.")
    print("The final expression is derived from the principles of magnetostatics for the given setup.")
    print("-" * 50)
    print("Final Expression for the change in mutual inductance per unit length (ΔM/L):")
    # Using unicode for subscripts and superscripts for better readability
    print("\u0394M/L = (\u03BC\u2080 * h\u00b2 * (R\u2082 - R\u2081)) / (2 * \u03c0 * d\u00b2 * R\u2081)")
    print("-" * 50)

    # --- Calculation with given numbers ---
    try:
        # Calculate the change in mutual inductance per unit length
        delta_M_per_L = (mu_0 * h**2 * (R2 - R1)) / (2 * math.pi * d**2 * R1)

        # --- Output the final equation with numbers ---
        print("Substituting the given values into the expression:")
        print(f"\u0394M/L = (({mu_0:.4e} H/m) * ({h} m)\u00b2 * ({R2} m - {R1} m)) / (2 * \u03c0 * ({d} m)\u00b2 * ({R1} m))")
        
        # Breaking down the calculation for clarity
        numerator = mu_0 * h**2 * (R2 - R1)
        denominator = 2 * math.pi * d**2 * R1
        print(f"\u0394M/L = {numerator:.4e} / {denominator:.4e}")
        
        print("\nResult:")
        print(f"\u0394M/L = {delta_M_per_L:.4e} H/m")

    except ZeroDivisionError:
        print("Error: R1 or d cannot be zero.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_inductance_change()