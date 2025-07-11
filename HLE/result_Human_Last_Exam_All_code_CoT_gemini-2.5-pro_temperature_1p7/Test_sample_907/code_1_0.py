def print_absorption_equation():
    """
    This script prints the equations for the absorption cross-section of a molecular chain
    under two different conditions: (a) no intermolecular interaction and (b) near-neighbor interaction.

    The equations are derived from first-order time-dependent perturbation theory for a Gaussian laser pulse.
    """
    
    # --- Define symbols for the equation ---
    sigma_omega = "σ(ω)"
    C = "C"  # Proportionality constant
    N = "N"  # Number of molecules in the chain
    mu_eg_sq = "|μ_eg|²"  # Squared transition dipole moment of a single molecule
    omega = "ω"  # Angular frequency of the laser
    omega_eg = "ω_eg"  # Transition frequency of a single molecule (E_eg / ħ)
    tau_sq = "τ²"  # Squared duration of the Gaussian pulse
    J_over_hbar = "J/ħ"  # Nearest-neighbor coupling energy in units of frequency

    # --- Case (a): No interaction between molecules ---
    print("Case a) The interaction between molecules can be neglected.")
    print("=" * 60)
    print("In this case, all N molecules absorb light independently. The total cross-section is the sum of the individual cross-sections. The absorption peak is centered at the single-molecule transition frequency ω_eg.")
    
    # Construct the equation string
    prefactor_a = f"{C} * {N} * {mu_eg_sq}"
    exponent_a = f"-({omega_eg} - {omega})² * {tau_sq}"
    equation_a = f"{sigma_omega} = {prefactor_a} * exp[{exponent_a}]"
    
    print("\nThe equation for the absorption cross-section is:")
    print(equation_a)
    print("\n")

    # --- Case (b): Near-neighbor interaction ---
    print("Case b) The interaction between near-neighbors is considered.")
    print("=" * 60)
    print("In this case, the molecular excitations are coupled, forming delocalized exciton states. Due to optical selection rules for a linear chain, only the k=0 exciton state is optically active. This shifts the absorption peak by an amount related to the coupling energy J and concentrates the oscillator strength into a single peak.")

    # Construct the equation string
    prefactor_b = f"{C} * {N} * {mu_eg_sq}"
    shifted_freq = f"({omega_eg} + 2*{J_over_hbar})"
    exponent_b = f"-({shifted_freq} - {omega})² * {tau_sq}"
    equation_b = f"{sigma_omega} = {prefactor_b} * exp[{exponent_b}]"

    print("\nThe equation for the absorption cross-section is:")
    print(equation_b)
    print("\nWhere:")
    print(f"  {sigma_omega:<10} = Absorption cross-section as a function of laser frequency ω")
    print(f"  {C:<10} = A constant of proportionality")
    print(f"  {N:<10} = The number of molecules in the chain")
    print(f"  {mu_eg_sq:<10} = The square of the transition dipole moment for a single molecule")
    print(f"  {omega_eg:<10} = The transition angular frequency of an isolated molecule")
    print(f"  {omega:<10} = The angular frequency of the incident laser pulse")
    print(f"  τ:<10} = The duration of the Gaussian laser pulse")
    print(f"  {J_over_hbar:<10} = The near-neighbor coupling energy (J) divided by ħ")


if __name__ == '__main__':
    print_absorption_equation()