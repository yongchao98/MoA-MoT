import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.

    The problem is solved using thin airfoil theory and the method of images.
    The final derived formula for the lift ratio is (1+K)/(1-K),
    where K is an interaction constant given by:
    K = (c/2) * (1/s - s / (s^2 + 4*h^2))

    Given values:
    - Separation s = 1/2 * c
    - Height h = 1/2 * c
    """

    # We can set the chord c=1.0 for simplicity, as it cancels out in the final ratio.
    c = 1.0

    # Define the separation (s) and ride height (h) based on the chord (c).
    s = 0.5 * c
    h = 0.5 * c

    print("Step 1: Define system geometry based on chord length c.")
    print(f"For c = {c:.1f}:")
    print(f"  - Separation, s = c/2 = {s:.2f}")
    print(f"  - Ride Height, h = c/2 = {h:.2f}")
    print("-" * 30)

    # Calculate the interaction constant K.
    print("Step 2: Calculate the interaction constant K.")
    print("Formula: K = (c/2) * (1/s - s/(s^2 + 4*h^2))")
    
    # Calculate intermediate terms for clarity.
    term1 = 1 / s
    term2 = s / (s**2 + 4 * h**2)
    K = (c / 2) * (term1 - term2)
    
    print(f"Value of K = ({c:.1f}/2) * (1/{s:.2f} - {s:.2f}/({s:.2f}^2 + 4*{h:.2f}^2))")
    print(f"Value of K = {c/2:.2f} * ({term1:.2f} - {term2:.2f}) = {K:.2f}")
    print("-" * 30)

    # Calculate the final lift ratio L1/L2.
    print("Step 3: Calculate the final lift ratio L1/L2.")
    print("Formula: L1/L2 = (1 + K) / (1 - K)")
    
    if abs(1 - K) < 1e-9:
        print("Error: Division by zero. The trailing aerofoil has zero lift.")
    else:
        lift_ratio = (1 + K) / (1 - K)
        numerator = 1 + K
        denominator = 1 - K
        print("\nFinal Equation with calculated values:")
        print(f"L1/L2 = (1 + {K:.2f}) / (1 - {K:.2f}) = {numerator:.2f} / {denominator:.2f} = {lift_ratio:.1f}")

# Execute the calculation function
calculate_lift_ratio()
<<<9.0>>>