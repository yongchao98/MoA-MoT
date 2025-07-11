def evaluate_quantum_correction():
    """
    This script provides a step-by-step derivation of the quantum correction
    to conductivity (weak localization) in a 3D bulk semiconductor.
    """
    print("Evaluating the quantum correction to conductivity in 3D (Weak Localization)\n")

    print("--- Step 1: The Formal Expression ---")
    print("The quantum correction to conductivity, Delta-sigma, arises from summing the 'cooperon' diagrams.")
    print("The expression in 3D is:")
    print("Delta-sigma = - (2 * e^2 * D / hbar) * Integral[d^3q / (2*pi)^3] * (1 / (D*q^2 + 1/tau_phi))\n")
    print("Where:")
    print("  e       = elementary charge")
    print("  hbar    = reduced Planck constant")
    print("  D       = electron diffusion constant")
    print("  tau_phi = phase coherence time")
    print("The integral over momentum 'q' is cut off at q_max = 1/l, where 'l' is the elastic mean free path.\n")

    print("--- Step 2: Simplification for Integration ---")
    print("We express the integral in spherical coordinates (d^3q -> 4*pi*q^2 dq) and simplify the prefactor:")
    print("Prefactor * Integral_angular_part = - (2*e^2*D/hbar) * (4*pi / (2*pi)^3) = - (e^2*D) / (pi^2 * hbar)")
    print("The integral becomes:")
    print("Delta-sigma = - e^2*D / (pi^2 * hbar) * Integral[from 0 to 1/l] dq * q^2 / (D*q^2 + 1/tau_phi)")
    print("\nTo simplify further, we introduce the phase coherence length, L_phi = sqrt(D * tau_phi):")
    print("Delta-sigma = - e^2 / (pi^2 * hbar) * Integral[from 0 to 1/l] dq * q^2 / (q^2 + 1/L_phi^2)\n")

    print("--- Step 3: Performing the Integration ---")
    print("Let a = 1/L_phi. We need to solve the integral of q^2 / (q^2 + a^2).")
    print("The standard form of this integral is: Integral[dq * q^2/(q^2+a^2)] = q - a * arctan(q/a).")
    print("Evaluating the definite integral from q=0 to q=1/l:")
    print("Result = [q - a * arctan(q/a)] from 0 to 1/l")
    print("Result = (1/l - (1/L_phi) * arctan((1/l) / (1/L_phi))) - (0 - a*arctan(0))")
    print("Result = 1/l - (1/L_phi) * arctan(L_phi / l)\n")

    print("--- Step 4: The Final Formula ---")
    print("Substituting the integral's result back into the expression for Delta-sigma:")
    print("Delta-sigma = - e^2 / (pi^2 * hbar) * [ 1/l - (1/L_phi) * arctan(L_phi / l) ]")
    print("This is often rearranged for clarity:")
    print("Delta-sigma = e^2 / (pi^2 * hbar) * [ (1/L_phi) * arctan(L_phi / l) - 1/l ]\n")

    print("--- Summary of Quantities in the Final Equation ---")
    print("e      : elementary charge")
    print("pi     : the mathematical constant ~3.14159")
    print("hbar   : reduced Planck constant (h / 2pi)")
    print("l      : elastic mean free path of the electron")
    print("L_phi  : phase coherence length of the electron")
    print("arctan : the inverse tangent function")

if __name__ == '__main__':
    evaluate_quantum_correction()