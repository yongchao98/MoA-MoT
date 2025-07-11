import sys

# This script explains the fundamental limit on the chemical potential for bosons
# in the context of Bose-Einstein condensation (BEC).

# The core of the explanation lies in the properties of the Bose-Einstein distribution function.
# We will print the step-by-step reasoning.

def explain_bec_chemical_potential():
    """
    Prints a step-by-step explanation for the limit on a boson's chemical potential.
    """
    print("Step 1: The Bose-Einstein distribution gives the average number of particles, <n_i>, in an energy state, ε_i:")
    print("  <n_i> = 1 / (exp((ε_i - μ) / (k_B * T)) - 1)\n")

    print("Step 2: For <n_i> to be a positive physical quantity, the denominator must be positive.")
    print("  This means: exp((ε_i - μ) / (k_B * T)) > 1\n")

    print("Step 3: Taking the natural log of both sides, we get the fundamental constraint on the chemical potential, μ:")
    print("  (ε_i - μ) > 0  which implies that  μ < ε_i\n")

    print("Step 4: This inequality must hold for all energy states. Therefore, μ must be less than the lowest possible energy state, the ground state energy, ε_0.")
    print("  This establishes the upper bound for the chemical potential: μ < ε_0\n")

    print("Step 5: For BEC to occur, a macroscopic number of particles must occupy the ground state. For the number of particles in the ground state, <n_0>, to be very large, the denominator must approach zero.")
    print("  <n_0> = 1 / (exp((ε_0 - μ) / (k_B * T)) - 1)")
    print("  This forces the chemical potential μ to get extremely close to ε_0 from below.\n")
    
    print("Step 6: At and below the critical temperature, μ becomes pinned at this limiting value. At absolute zero (T=0), all particles are in the ground state, and the chemical potential is exactly the ground state energy.")
    print("Final Limit Equation:")
    print("  μ = ε_0")
    print("  where ε_0 is the ground state energy, which is also equal to the chemical potential of a non-interacting Bose gas at zero temperature.\n")

explain_bec_chemical_potential()
sys.stdout.flush()