import sys

def solve_absorption_cross_section():
    """
    This script explains and provides the equations for the absorption cross-section
    of a molecular chain interacting with an ultrashort Gaussian laser pulse,
    based on first-order time-dependent perturbation theory.
    """

    # --- Introduction to the General Theory ---
    print("--- General Theoretical Framework ---")
    print("The absorption cross-section σ(ω_L) for a system absorbing a laser pulse with carrier frequency ω_L")
    print("can be derived using first-order time-dependent perturbation theory. The general formula is")
    print("proportional to the sum over all possible final states |f> of the transition probability from")
    print("the initial state |i> to the final state |f>.")
    print("\nFor an ultrashort Gaussian-shaped pulse of duration τ, the theory gives the following expression:")
    print("\n  σ(ω_L) ∝ Σ_f  ω_fi * |<f|μ|i>|^2 * exp(-(ω_fi - ω_L)^2 * τ^2)\n")
    print("where:")
    print("  - ω_L is the laser carrier frequency.")
    print("  - |i> is the initial (ground) state of the molecular chain.")
    print("  - |f> are the final (excited) states.")
    print("  - ħω_fi = (E_f - E_i) is the transition energy from state |i> to |f>.")
    print("  - μ is the electric dipole moment operator of the entire chain.")
    print("  - τ is related to the duration of the Gaussian laser pulse.")
    print("  - ħ is the reduced Planck constant.")
    print("\nLet's apply this general formula to the two specific cases.")
    print("-" * 50)

    # --- Case a) No Interaction ---
    print("\na) Case: The interaction between molecules can be neglected.\n")
    print("In this scenario, the N molecules in the chain are treated as independent entities.")
    print("The final excited states are N-fold degenerate, where each state corresponds to a single")
    print("molecule being excited from its ground state. All these transitions have the same properties:")
    print("  - The transition frequency is ω_0, the natural frequency of a single molecule.")
    print("  - The transition dipole moment is μ_0, the dipole moment of a single molecule.")
    print("\nThe total cross-section is the incoherent sum of the contributions from each of the N molecules.")
    print("The sum Σ_f |<f|μ|i>|^2 becomes N * |μ_0|^2.")
    print("\nThe resulting equation for the absorption cross-section is:")
    print("\n  σ(ω_L) = C * N * ω_0 * |μ_0|^2 * exp(-(ω_0 - ω_L)^2 * τ^2)\n")
    print("where:")
    print("  - C is a constant of proportionality.")
    print("  - N is the number of molecules in the chain.")
    print("  - The spectrum consists of a single Gaussian peak at the monomer frequency ω_0.")
    print("-" * 50)

    # --- Case b) Nearest-Neighbor Interaction ---
    print("\nb) Case: The interaction between near-neighbors should be considered.\n")
    print("Here, we use the Frenkel exciton model. The coupling J between nearest neighbors mixes the")
    print("excited states of individual molecules into collective, delocalized 'exciton' states.")
    print("These exciton states have different energies, forming an 'exciton band'.")
    print("\nFor a chain with periodic boundary conditions, a crucial optical selection rule emerges:")
    print("only the exciton state with wavevector k=0 is optically active (bright).")
    print("  - The transition frequency to this state is shifted by the interaction: ω_k=0 = ω_0 + 2*J/ħ.")
    print("  - The transition dipole moment is coherently enhanced: |μ_k=0|^2 = N * |μ_0|^2.")
    print("\nSince only one transition is allowed, the sum Σ_f collapses to a single term.")
    print("The resulting equation for the absorption cross-section is:")
    print("\n  σ(ω_L) = C * N * (ω_0 + 2*J/ħ) * |μ_0|^2 * exp(-(ω_0 + 2*J/ħ - ω_L)^2 * τ^2)\n")
    print("where:")
    print("  - J is the nearest-neighbor coupling energy (can be positive or negative).")
    print("  - The spectrum is a single Gaussian peak, but its position is shifted by 2*J/ħ")
    print("    relative to the non-interacting case. This is the exciton shift.")
    print("-" * 50)

if __name__ == '__main__':
    solve_absorption_cross_section()