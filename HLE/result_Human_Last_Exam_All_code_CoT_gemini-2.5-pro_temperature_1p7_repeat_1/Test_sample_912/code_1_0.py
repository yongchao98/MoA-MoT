import sympy as sp

def solve_magnetic_work_problem():
    """
    This function provides a step-by-step derivation for the work done
    by the current source in the given magnetic circuit problem.
    """

    # --- Introduction ---
    print("This script will derive the formula for the work done by the current source in one cycle.")
    print("The work is calculated by integrating the electrical power over the cycle.")
    print("-" * 50)

    # --- Step 1: System Inductance L(x) ---
    print("Step 1: Determine the inductance L(x) as a function of the block's position x.")
    print("The magnetic yoke has infinite permeability, so its reluctance is zero.")
    print("The total reluctance of the circuit is determined by the air gap. The gap consists of two parallel magnetic paths:")
    print("  - Path 1: Through the magnetic block, with area w*x and permeability mu.")
    print("    Reluctance R_block = g / (mu * w * x)")
    print("  - Path 2: Through the remaining air, with area w*(D-x) and permeability mu_0.")
    print("    Reluctance R_air = g / (mu_0 * w * (D - x))")
    print("The total reluctance R(x) is the parallel combination of R_block and R_air:")
    print("  R(x) = (1/R_block + 1/R_air)^-1 = g / (w * ((mu - mu_0)*x + mu_0*D))")
    print("The inductance L(x) is given by N^2 / R(x):")
    print("  L(x) = (N^2 * w / g) * ((mu - mu_0)*x + mu_0*D)")
    print("-" * 50)

    # --- Step 2: Work Calculation over the Cycle ---
    print("Step 2: Calculate the work done by the source over the rectangular cycle.")
    print("The incremental work done by the source is dW = I * d(lambda), where lambda = L(x)*I is the flux linkage.")
    print("We integrate this over the four steps of the cycle (A->B, B->C, C->D, D->A).\n")
    print("  A=(x1, I1), B=(x2, I1), C=(x2, I2), D=(x1, I2)\n")

    print("  Step A->B (x moves from x1 to x2 at constant I=I1):")
    print("    W_AB = integral from x1 to x2 of I1^2 * dL/dx dx = I1^2 * (L(x2) - L(x1))")

    print("  Step B->C (I increases from I1 to I2 at constant x=x2):")
    print("    W_BC = integral from I1 to I2 of L(x2) * I dI = (1/2) * L(x2) * (I2^2 - I1^2)")

    print("  Step C->D (x moves from x2 to x1 at constant I=I2):")
    print("    W_CD = integral from x2 to x1 of I2^2 * dL/dx dx = I2^2 * (L(x1) - L(x2))")

    print("  Step D->A (I decreases from I2 to I1 at constant x=x1):")
    print("    W_DA = integral from I2 to I1 of L(x1) * I dI = (1/2) * L(x1) * (I1^2 - I2^2)")
    print("-" * 50)

    # --- Step 3: Summing the Work ---
    print("Step 3: Sum the work from all four steps.")
    print("  W_cycle = W_AB + W_BC + W_CD + W_DA")
    print("  W_cycle = I1^2*(L(x2)-L(x1)) + 0.5*L(x2)*(I2^2-I1^2) - I2^2*(L(x2)-L(x1)) - 0.5*L(x1)*(I2^2-I1^2)")
    print("  Combining terms:")
    print("  W_cycle = (L(x2) - L(x1)) * (I1^2 - I2^2) + 0.5 * (L(x2) - L(x1)) * (I2^2 - I1^2)")
    print("  W_cycle = -(L(x2) - L(x1)) * (I2^2 - I1^2) + 0.5 * (L(x2) - L(x1)) * (I2^2 - I1^2)")
    print("  W_cycle = -0.5 * (L(x2) - L(x1)) * (I2^2 - I1^2)")
    print("-" * 50)

    # --- Step 4: Final Expression ---
    print("Step 4: Substitute the expression for L(x) into the result.")
    print("  From Step 1, L(x2) - L(x1) = (N^2 * w * (mu - mu_0) / g) * (x2 - x1)")
    print("  Substituting this into the equation for W_cycle gives the final answer.")

    print("\n--- Final Equation ---")
    print("The work W done by the current source for each cycle is:")

    # Printing each "number" or component of the final equation as requested.
    print("W = ", end="")
    print("-(1/2) ", end="")
    print("* [ (mu - mu_0) * N^2 * w / g ] ", end="")
    print("* (x_2 - x_1) ", end="")
    print("* (I_2^2 - I_1^2)")
    
    print("\nThis can be written as:")
    print("W = - ( (mu - mu_0) / (2 * g) ) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)")

if __name__ == '__main__':
    solve_magnetic_work_problem()