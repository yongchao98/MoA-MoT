import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem near the ground.
    
    The problem is solved using thin aerofoil theory and the mirror image method.
    The lift ratio is found to be (A+C)/(A-C), where A and C are coefficients
    derived from the geometry and vortex interaction model.
    """
    
    # We can use non-dimensional quantities. Let chord c = 1.0.
    c = 1.0
    
    # Given parameters:
    # Separation distance s = 1/2 * c
    s = 0.5 * c
    # Ride height h = c / 2
    h = 0.5 * c
    
    print("Problem Parameters (assuming chord c=1):")
    print(f"Separation distance, s = {s}")
    print(f"Ride height, h = {h}\n")
    
    # Step 1: Calculate the coefficient A
    # A = 1 - c / (4 * h)
    A = 1 - c / (4 * h)
    print("Calculating coefficient A:")
    print(f"A = 1 - c / (4 * h) = 1 - {c} / (4 * {h}) = {A}\n")
    
    # Step 2: Calculate the coefficient C
    # C = (c / 2) * (1/s - s / (s**2 + 4 * h**2))
    term1_C = 1 / s
    term2_C = s / (s**2 + 4 * h**2)
    C = (c / 2) * (term1_C - term2_C)
    print("Calculating coefficient C:")
    print(f"C = (c / 2) * [1/s - s / (s^2 + 4*h^2)]")
    print(f"C = ({c} / 2) * [1/{s} - {s} / ({s}^2 + 4*{h}^2)]")
    print(f"C = 0.5 * [{term1_C} - {s} / ({s**2} + {4*h**2})]")
    print(f"C = 0.5 * [{term1_C} - {term2_C}] = {C}\n")

    # Step 3: Calculate the lift ratio L1/L2
    # L1/L2 = (A + C) / (A - C)
    if (A - C) == 0:
        lift_ratio = float('inf')
    else:
        lift_ratio = (A + C) / (A - C)

    print("Final Calculation for the Lift Ratio L1/L2:")
    print(f"L1/L2 = (A + C) / (A - C)")
    print(f"L1/L2 = ({A} + {C}) / ({A} - {C})")
    print(f"L1/L2 = {A + C} / {A - C}")
    print(f"L1/L2 = {lift_ratio}")

calculate_lift_ratio()