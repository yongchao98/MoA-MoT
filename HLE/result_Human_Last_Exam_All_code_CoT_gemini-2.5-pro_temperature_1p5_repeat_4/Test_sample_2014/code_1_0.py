import numpy as np

def solve_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem formation
    in ground effect using the mirror image method.
    """

    # Step 1: Define Problem Parameters
    # We are given s = 1/2 * c and h = 1/2 * c.
    # We normalize the chord length c to 1.0 for simplicity.
    c = 1.0
    h = 0.5 * c
    s = 0.5 * c

    print("--- Problem Setup ---")
    print(f"Normalized Chord Length c = {c:.1f}")
    print(f"Ride Height h = {h:.1f}")
    print(f"Tandem Separation s = {s:.1f}\n")

    # Step 2: Formulate the System of Equations
    # From thin aerofoil theory and the vortex lattice method, we can derive a
    # system of linear equations for the circulations (Gamma) or lifts (L)
    # of the two aerofoils. The lift L is directly proportional to circulation.
    #
    # Assuming the same geometric angle of attack for both aerofoils, we get:
    # K = A * L1 + B * L2
    # K = C * L1 + D * L2
    # where K is a constant related to the angle of attack and freestream velocity.
    # The coefficients A, B, C, D depend on the geometry (h, s, c).

    # The analytical derivation yields the following coefficients:
    # A = (1/c + C_11) = 1 + 1/5 = 6/5
    # B = C_12 = 1/4
    # C = C_21 = -1/20
    # D = (1/c + C_22) = 1 + 1/5 = 6/5

    A = 1.2
    B = 0.25
    C = -0.05
    D = 1.2

    print("--- System of Equations for Lifts L1 and L2 ---")
    print("Let K be a constant proportional to the angle of attack.")
    print(f"Equation for Aerofoil 1: K = ({A:.3f}) * L1 + ({B:.3f}) * L2")
    print(f"Equation for Aerofoil 2: K = ({C:.3f}) * L1 + ({D:.3f}) * L2\n")

    # Step 3: Solve for the Lift Ratio L1/L2
    # By equating the two expressions for K, we can solve for the ratio L1/L2.
    # A * L1 + B * L2 = C * L1 + D * L2
    # (A - C) * L1 = (D - B) * L2
    # L1 / L2 = (D - B) / (A - C)

    coeff_L1 = A - C
    coeff_L2 = D - B

    lift_ratio = coeff_L2 / coeff_L1

    print("--- Solving for the Lift Ratio ---")
    print("Equating the two expressions gives:")
    print(f"({A:.3f})*L1 + ({B:.3f})*L2 = ({C:.3f})*L1 + ({D:.3f})*L2")
    print("Rearranging the terms yields:")
    print(f"({A:.3f} - ({C:.3f})) * L1 = ({D:.3f} - {B:.3f}) * L2")
    print(f"({coeff_L1:.3f}) * L1 = ({coeff_L2:.3f}) * L2")
    print("\nThe final lift ratio is:")
    print(f"L1 / L2 = {coeff_L2:.3f} / {coeff_L1:.3f} = {lift_ratio:.4f}")

    return lift_ratio

if __name__ == '__main__':
    solve_lift_ratio()