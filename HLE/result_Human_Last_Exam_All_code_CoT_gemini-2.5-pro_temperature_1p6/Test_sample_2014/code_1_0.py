def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the mirror image method.
    """
    # --- 1. Explanation of the Method ---
    print("This solution models the aerofoils as vortices and uses the mirror image method for ground effect.")
    print("\n--- Model and Equations ---")
    print("Aerofoil 1 (front): Circulation Γ1 at position (0, h)")
    print("Aerofoil 2 (rear):  Circulation Γ2 at position (s, h)")
    print("Image Vortices: -Γ1 at (0, -h) and -Γ2 at (s, -h)")
    print("\nGiven parameters:")
    print("Aerofoil chord = c")
    print("Separation, s = c / 2")
    print("Ride height, h = c / 2")

    print("\nThe flow tangency condition leads to a system of equations for identical aerofoils:")
    print("  1) Γ1 + A * Γ2 = K")
    print("  2) Γ2 - A * Γ1 = K")
    print("Where 'A' is an interaction factor defined as: A = (c/2) * [1/s - s/(s² + 4h²)]")

    print("\nSolving this system for the lift ratio (L1/L2 = Γ1/Γ2) yields:")
    print("  L1/L2 = (1 - A) / (1 + A)")

    print("\n--- Numerical Calculation ---")
    # We can use dimensionless values by setting c=1, which makes s=0.5 and h=0.5
    c = 1.0
    s = 0.5 * c
    h = 0.5 * c
    print(f"To calculate, let c = {c:.1f}. This means s = {s:.1f} and h = {h:.1f}.")

    # Calculate the interaction factor A
    A_val = (c / 2.0) * (1.0/s - s / (s**2 + 4 * h**2))
    
    print("\nFirst, we calculate the value of A:")
    print(f"A = ({c:.1f} / 2) * [1/{s:.1f} - {s:.1f} / ({s:.1f}² + 4*{h:.1f}²)]")
    print(f"A = 0.5 * [2.0 - {s:.1f} / ({s**2:.2f} + {4 * h**2:.2f})]")
    print(f"A = 0.5 * [2.0 - {s:.1f} / {s**2 + 4 * h**2:.2f}]")
    print(f"A = 0.5 * [2.0 - {s / (s**2 + 4 * h**2):.1f}]")
    print(f"A = 0.5 * 1.6")
    print(f"A = {A_val:.1f}")

    # Calculate the final lift ratio L1/L2
    numerator = 1.0 - A_val
    denominator = 1.0 + A_val
    ratio = numerator / denominator

    print("\nFinally, we calculate the lift ratio L1/L2 using the value of A:")
    print(f"L1/L2 = (1 - A) / (1 + A)")
    print(f"L1/L2 = (1 - {A_val:.1f}) / (1 + {A_val:.1f})")
    print(f"L1/L2 = {numerator:.1f} / {denominator:.1f}")
    
    print("\n--------------------")
    print(f"The final lift ratio is: {ratio:.4f}")
    print("As a fraction, the result is 0.2 / 1.8 = 2 / 18 = 1/9.")
    print("--------------------")

# Execute the calculation function
calculate_lift_ratio()
