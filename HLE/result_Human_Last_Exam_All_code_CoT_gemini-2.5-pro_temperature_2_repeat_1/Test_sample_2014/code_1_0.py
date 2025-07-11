import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem and in ground effect.

    The problem parameters are:
    - Separation s = c/2
    - Ride height h = c/2

    The function derives and calculates the lift ratio based on potential flow theory
    and the mirror image method.
    """

    # We can set the chord 'c' to any value since it will cancel out.
    # Let's choose c = 2.0 for simpler intermediate numbers.
    c = 2.0
    s = 0.5 * c
    h = 0.5 * c

    # --- Explanation ---
    print("This script calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.")
    print("The model uses point vortices and the mirror image method.")
    print("The key formulas are derived from the flow tangency conditions for each aerofoil.")
    print("\nDerivation leads to the lift ratio formula:")
    print("L1/L2 = (2/c + K') / (2/c - K')")
    print("where K' is the interaction factor: K' = 1/s - s / (s^2 + 4h^2)")
    print("-" * 30)

    # --- Step-by-step Calculation ---
    print("\n1. Given Parameters:")
    print(f"   Chord c")
    print(f"   Separation s = 1/2 * c")
    print(f"   Ride height h = c/2")

    print("\n2. Calculate the interaction factor K' in terms of c:")
    # K' = 1/(c/2) - (c/2) / ((c/2)^2 + 4*(c/2)^2)
    #    = 2/c - (c/2) / (c^2/4 + 4*c^2/4)
    #    = 2/c - (c/2) / (5c^2/4)
    #    = 2/c - (c/2) * (4/5c^2)
    #    = 2/c - 2/(5c)
    #    = (10 - 2) / (5c) = 8/(5c)
    print("   Substituting s=c/2 and h=c/2 into the formula for K':")
    print("   K' = 1/(c/2) - (c/2) / ((c/2)^2 + 4*(c/2)^2)")
    print("   K' = 2/c - (c/2) / (5c^2/4)")
    print("   K' = 2/c - 2/(5c)")
    print("   K' = 8/(5c)")

    print("\n3. Substitute K' into the lift ratio formula:")
    # Using the symbolic result for K'
    # L1/L2 = (2/c + 8/(5c)) / (2/c - 8/(5c))
    # We can factor out 1/c from the numerator and denominator
    # L1/L2 = (2 + 8/5) / (2 - 8/5)
    #       = (10/5 + 8/5) / (10/5 - 8/5)
    #       = (18/5) / (2/5)
    #       = 18 / 2 = 9
    print("   L1/L2 = (2/c + 8/(5*c)) / (2/c - 8/(5*c))")
    term1_num = 2
    term2_num = "8/5"
    final_num = "18/5"
    
    term1_den = 2
    term2_den = "8/5"
    final_den = "2/5"

    ratio_numerator = 18 / 5
    ratio_denominator = 2 / 5
    final_ratio = ratio_numerator / ratio_denominator
    
    print(f"   L1/L2 = ( {term1_num} + {term2_num} ) / ( {term1_den} - {term2_den} )   (after canceling '1/c')")
    print(f"   L1/L2 = ( {final_num} ) / ( {final_den} )")
    print(f"   L1/L2 = {ratio_numerator} / {ratio_denominator}")

    print("\n4. Final Result:")
    print(f"   L1/L2 = {final_ratio}")
    print("-" * 30)

# Run the calculation and print the results
calculate_lift_ratio()
print("<<<9.0>>>")