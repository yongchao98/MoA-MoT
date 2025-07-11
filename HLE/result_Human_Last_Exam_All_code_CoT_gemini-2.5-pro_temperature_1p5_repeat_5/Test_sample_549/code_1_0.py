import math

def evaluate_quantum_correction_3d():
    """
    This script provides a step-by-step evaluation of the quantum correction 
    to conductivity (weak localization) in a 3D bulk semiconductor.
    """
    
    # 1. Explain the setup.
    print("### Evaluation of Quantum Correction to Conductivity in 3D ###")
    print("This evaluation derives the formula for the weak localization correction, which is a quantum interference effect.")
    print("The correction is proportional to the probability of an electron returning to its starting point during diffusion.\n")

    # 2. Start the derivation.
    print("--- Step 1: Probability of a Diffusing Particle in 3D ---")
    print("The probability density P(r, t) for a particle at position r at time t, starting from the origin, is given by the diffusion equation solution:")
    print("  P(r, t) = [1 / (4 * pi * D * t)^(3/2)] * exp[-r^2 / (4 * D * t)]")
    print("where D is the diffusion coefficient.\n")

    print("--- Step 2: Probability of Returning to the Origin (r=0) ---")
    print("To find the probability density of the electron returning to its starting point, we set r=0:")
    print("  P(0, t) = [1 / (4 * pi * D * t)^(3/2)]\n")

    print("--- Step 3: Integrate Over Relevant Timescales ---")
    print("We integrate this probability density over the time during which quantum interference is possible.")
    print("  - Lower limit: The momentum relaxation time (tau), the shortest time for diffusion to be valid.")
    print("  - Upper limit: The phase coherence time (tau_phi), after which interference is destroyed.")
    print("The total integrated return probability density (W) is:")
    print("  W = Integral[from tau to tau_phi] P(0, t) dt\n")

    print("--- Step 4: Change of Variables for Simplification ---")
    print("To simplify the integral, we switch from time 't' to diffusion length 'L', using the relationship L^2 = D*t.")
    print("This gives: t = L^2 / D, and the differential dt = (2*L / D) dL.")
    print("The integration limits transform to:")
    print("  - Lower limit: Elastic mean free path, l = sqrt(D * tau)")
    print("  - Upper limit: Phase coherence length, L_phi = sqrt(D * tau_phi)")
    print("The term (D*t)^(3/2) in P(0,t) becomes L^3. The integral becomes a path integral:\n")

    print("--- Step 5: Solve the Transformed Integral ---")
    print("  W = Integral[from l to L_phi] (1 / (4*pi*L^3)) * (8*pi*L/D) dL ... No, the proportionality is simpler.")
    print("The return probability integral is proportional to Integral[ (1/L^3) * (L/D) dL ]")
    print("  W_proportional = Integral[from l to L_phi] L^(-2) dL")
    print("Solving this elementary integral:")
    print("  W_proportional = [-1/L] from l to L_phi")
    print("                 = [(-1/L_phi) - (-1/l)]")
    print("                 = (1/l - 1/L_phi)\n")
    
    print("--- Step 6: Relate Return Probability to Conductivity Correction ---")
    print("A full diagrammatic calculation (summing Cooperon loops) shows that the conductivity correction, delta_sigma, is:")
    print("  delta_sigma = - C_3D * (e^2 / hbar) * (1/l - 1/L_phi)")
    print("where 'e' is the electron charge, 'hbar' is the reduced Planck constant, and C_3D is a dimensionless prefactor from the full theory.\n")

    print("--- Step 7: Final Expression for the Quantum Correction ---")
    print("From the rigorous theory of weak localization, the prefactor C_3D for a 3D system is found to be 1 / (2 * pi^2).")
    print("Therefore, the final result for the quantum correction to conductivity is:\n")
    
    final_formula = "delta_sigma = - (e^2 / (2 * pi^2 * hbar)) * (1/l - 1/L_phi)"
    print("=====================================================================")
    print(f"  {final_formula}")
    print("=====================================================================\n")

    print("Here is a breakdown of the numbers and symbols in the final equation:")
    print("  - The equation is NEGATIVE, as weak localization reduces conductivity.")
    print("  - 'e' represents the elementary charge.")
    print("  - The number '2' is in the denominator.")
    print("  - 'pi' is the mathematical constant (approx 3.14159), which appears squared.")
    print("  - 'hbar' is the reduced Planck constant.")
    print("  - 'l' is the elastic mean free path.")
    print("  - 'L_phi' is the phase coherence length.")

if __name__ == '__main__':
    evaluate_quantum_correction_3d()