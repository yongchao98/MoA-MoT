def solve_perturbation_theory():
    """
    Calculates the energy spectrum for a harmonic oscillator with a quartic perturbation,
    based on the provided self-energy diagram. The code prints the step-by-step derivation.
    """

    print("The energy spectrum of the perturbed system is found from the poles of the full propagator D(ω).")
    print("The full propagator is related to the free propagator D₀(ω) and the self-energy Σ*(ω) by the Dyson equation:")
    print("D(ω) = D₀(ω) / (1 + i * Σ*(ω) * D₀(ω))")
    print("The poles are found where the denominator is zero, which leads to a condition for the new frequency ω'.\n")

    print("Step 1: Calculate the Self-Energy Σ*")
    print("---------------------------------------")
    print("The self-energy diagram corresponds to a vertex with a loop.")
    print("The interaction Hamiltonian is H_int = (u/4!) * x⁴.")
    print("The vertex rule combined with the combinatorial factor for this diagram gives a factor of u/2.")
    print("The self-energy Σ* is given by: Σ* = (u/2) * D₀(0)")
    print("where D₀(0) = <0|x(0)x(0)|0> is the free propagator at t=0.")
    print("For the ground state of the simple harmonic oscillator, this expectation value is:")
    print("D₀(0) = <0|x²|0> = ħ / (2 * m * ω₀)")
    print("\nSubstituting this in, we get:")
    print("Σ* = (u/2) * (ħ / (2 * m * ω₀))")
    print("Σ* = (u * ħ) / (4 * m * ω₀)\n")

    print("Step 2: Find the New Frequency ω'")
    print("---------------------------------")
    print("The pole of the full propagator D(ω) defines the new frequency ω'. The pole condition is:")
    print("D₀(ω')⁻¹ + i * Σ* = 0")
    print("The inverse free propagator is D₀(ω)⁻¹ = -i * (m/ħ) * (ω² - ω₀²).")
    print("So, -i * (m/ħ) * (ω'² - ω₀²) + i * Σ* = 0")
    print("This simplifies to: ω'² - ω₀² = (ħ/m) * Σ*")
    print("\nPlugging in our result for Σ*:")
    print("ω'² = ω₀² + (ħ/m) * [ (u * ħ) / (4 * m * ω₀) ]")
    print("ω'² = ω₀² + (u * ħ²) / (4 * m² * ω₀)\n")

    print("Step 3: Determine the Energy Spectrum Eₙ")
    print("-------------------------------------------")
    print("This self-energy calculation renormalizes the oscillator's frequency from ω₀ to ω'.")
    print("The system is now described as a new harmonic oscillator with frequency ω'.")
    print("The energy levels of a harmonic oscillator with frequency ω' are given by:")
    print("Eₙ = ħ * ω' * (n + 1/2)")
    print("\nTherefore, the predicted energy spectrum is:")
    print("\nEₙ = ħ * sqrt( ω₀² + (u * ħ²) / (4 * m² * ω₀) ) * ( n + 1/2 )\n")
    print("Final Equation Components:")
    print("Term 1 (Planck's constant): ħ")
    print("Term 2 (New frequency ω'): sqrt( ω₀² + (u * ħ²) / (4 * m² * ω₀) )")
    print("Term 3 (Quantum number part): (n + 1/2)")

solve_perturbation_theory()