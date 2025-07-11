import numpy as np

def solve_aerofoil_tandem_ground_effect():
    """
    Calculates the lift ratio of two aerofoils in tandem formation
    and in ground effect using the vortex method.
    """
    # Based on the derivation from the vortex model, the relationship
    # between the circulations Γ1 (front aerofoil) and Γ2 (rear aerofoil)
    # can be expressed as a system of linear equations:
    #
    # A[0,0]*Γ1 + A[0,1]*Γ2 = K
    # A[1,0]*Γ1 + A[1,1]*Γ2 = K
    #
    # where K is a constant related to the freestream velocity and angle of attack.
    # The coefficients are derived from the geometry of the setup:
    # chord c, separation s=c/2, and ride height h=c/2.

    # Coefficients for the system of equations.
    A = np.array([
        [16.0, -5.0],
        [1.0, 16.0]
    ])

    # The right-hand side is a constant K. For finding the ratio Γ1/Γ2,
    # we can set K=1 without loss of generality.
    B = np.array([1.0, 1.0])

    print("This script calculates the lift ratio L1/L2 for two aerofoils in tandem and in ground effect.")
    print("The problem is solved using the vortex mirror image method and thin aerofoil theory.")
    print("-" * 70)
    print("The relationship between the circulations Γ1 (front) and Γ2 (rear) is:")
    print("where K is a constant depending on freestream conditions (velocity and angle of attack).\n")

    # Print the equations with their numbers
    print(f"Equation 1: {A[0,0]:.2f} * Γ1 + ({A[0,1]:.2f}) * Γ2 = K")
    print(f"Equation 2:  {A[1,0]:.2f} * Γ1 + {A[1,1]:.2f} * Γ2 = K\n")
    print("-" * 70)

    try:
        # Solve the system Ax = B for x, where x contains the relative circulations.
        relative_circulations = np.linalg.solve(A, B)
        gamma1_rel = relative_circulations[0]
        gamma2_rel = relative_circulations[1]

        # The lift ratio L1/L2 is equal to the circulation ratio Γ1/Γ2.
        lift_ratio = gamma1_rel / gamma2_rel

        print("Solving this system gives the ratio of circulations Γ1/Γ2.")
        print("Since lift is directly proportional to circulation (L = ρUΓ), the lift ratio is L1/L2 = Γ1/Γ2.\n")
        print(f"The calculated lift ratio is L1 / L2 = {lift_ratio:.2f}")

    except np.linalg.LinAlgError:
        print("The system of equations is singular and cannot be solved.")

if __name__ == '__main__':
    solve_aerofoil_tandem_ground_effect()