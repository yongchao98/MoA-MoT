import math

def solve_sphere_radiation_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for the described rotating sphere.

    The final expression for the ratio is derived as:
    Ratio = (π/3) * (16 - 7√2)
    """

    # Define the constants from the derived formula
    pi = math.pi
    sqrt2 = math.sqrt(2)
    
    val_A = pi
    val_B = 3.0
    val_C = 16.0
    val_D = 7.0
    val_E = sqrt2
    
    # Calculate the final numerical result
    result = (val_A / val_B) * (val_C - val_D * val_E)

    # Output the explanation and the final answer
    print("The maximum achievable ratio is derived from the following equation:")
    print("Ratio = (A / B) * (C - D * E)\n")
    print("Where the component values are:")
    print(f"  A = π ≈ {val_A:.6f}")
    print(f"  B = {val_B}")
    print(f"  C = {val_C}")
    print(f"  D = {val_D}")
    print(f"  E = √2 ≈ {val_E:.6f}\n")
    
    # Show the final calculation
    term1 = val_A / val_B
    term2 = val_C - val_D * val_E
    print("Substituting these values into the equation:")
    print(f"Ratio = ({val_A:.6f} / {val_B}) * ({val_C} - {val_D} * {val_E:.6f})")
    print(f"Ratio = {term1:.6f} * {term2:.6f}")
    print(f"\nMaximum Ratio ≈ {result:.6f}")


solve_sphere_radiation_ratio()
<<<6.388533>>>