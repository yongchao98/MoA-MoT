import numpy as np

def solve_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem near the ground.
    
    This function follows these steps:
    1.  Sets up the problem using vortex modeling for the aerofoils and their images.
    2.  Derives the influence coefficients that form a system of linear equations for the vortex strengths (Gamma1, Gamma2).
    3.  Solves the system algebraically to find the ratio Gamma1/Gamma2.
    4.  Prints the final lift ratio, as L1/L2 = Gamma1/Gamma2.
    """

    print("--- Aerofoil Lift Ratio Calculation ---")
    print("This script calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.")
    print("Problem parameters:")
    print("  - Chord length: c")
    print("  - Separation (TE-LE): s = 0.5 * c")
    print("  - Ride Height: h = 0.5 * c\n")

    print("Step 1: Formulating the flow tangency equations.")
    print("The problem is modeled by a system of two linear equations:")
    print("  U*alpha = A * Gamma1 + B * Gamma2  (for Aerofoil 1)")
    print("  U*alpha = C * Gamma1 + D * Gamma2  (for Aerofoil 2)")
    print("where A, B, C, D are influence coefficients based on the geometry.\n")

    print("Step 2: Calculating the influence coefficients.")
    print("After substituting s=c/2 and h=c/2 into the downwash equations, the coefficients are:")
    print("  A = 1 / (5*pi*c)")
    print("  B = -3 / (4*pi*c)")
    print("  C = -1 / (20*pi*c)")
    print("  D = 1 / (5*pi*c)\n")

    print("Step 3: Setting up the system of equations.")
    print("Since U*alpha is the same for both, we can equate the right-hand sides:")
    print("  A*Gamma1 + B*Gamma2 = C*Gamma1 + D*Gamma2")
    print("Substituting the coefficients (and canceling pi*c):")
    print("  (1/5)*Gamma1 - (3/4)*Gamma2 = (-1/20)*Gamma1 + (1/5)*Gamma2\n")

    print("Step 4: Solving for the ratio Gamma1 / Gamma2.")
    print("Rearranging the terms to group Gamma1 and Gamma2:")
    print("  (1/5 + 1/20)*Gamma1 = (1/5 + 3/4)*Gamma2")
    # To work with integers, we can multiply the whole equation by 20.
    # (4+1)*Gamma1 = (4+15)*Gamma2
    g1_coeff_num = 4 + 1
    g2_coeff_num = 4 + 15
    print(f"  ({g1_coeff_num}/20)*Gamma1 = ({g2_coeff_num}/20)*Gamma2")
    print(f"  {g1_coeff_num}*Gamma1 = {g2_coeff_num}*Gamma2\n")

    print("Step 5: Final Result.")
    print("The lift ratio L1/L2 is equal to the circulation ratio Gamma1/Gamma2.")
    # From 5*Gamma1 = 19*Gamma2, we get Gamma1/Gamma2 = 19/5
    numerator = g2_coeff_num
    denominator = g1_coeff_num
    result = numerator / denominator
    
    print("The final equation is L1 / L2 = Gamma1 / Gamma2")
    print(f"L1 / L2 = {numerator} / {denominator}")
    print(f"L1 / L2 = {result}")

# Execute the function to get the answer
solve_lift_ratio()
<<<3.8>>>