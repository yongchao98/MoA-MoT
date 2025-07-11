import scipy.constants as const

def calculate_energy_difference():
    """
    This function solves for the energy difference ΔE based on the
    particle-on-a-sphere model derived from the problem's potential.
    """
    # Given parameters
    R_nm = 3.0  # nm
    m_particle = 9.11e-31  # kg

    # Convert parameters to SI units
    R = R_nm * 1e-9  # meters

    # Physical constants from scipy
    hbar = const.hbar      # Reduced Planck constant in J·s
    eV_to_J = const.e      # elementary charge in C, for J to eV conversion

    # As derived in the plan, the energy levels for a particle on a sphere
    # are E_l = l(l+1)ħ² / (2mR²).
    # E₁ is the ground state (l=0), E₂ is the first excited state (l=1).
    E1_J = 0
    E2_J = (hbar**2) / (m_particle * R**2)

    # The energy difference is ΔE = E₂ - E₁.
    delta_E_J = E2_J - E1_J
    delta_E_eV = delta_E_J / eV_to_J

    print("The potential suggests the particle is confined to a spherical shell of radius R.")
    print("This is modeled as a particle on a sphere, with energy levels E_l = l(l+1)ħ²/(2mR²).")
    print("\nFirst level E₁ (l=0) = 0 eV")
    print(f"Second level E₂ (l=1) = ħ² / (m * R²)")
    print(f"\nEnergy Difference ΔE = E₂ - E₁ = ħ² / (m * R²)")
    print("\nCalculation with numerical values:")
    # The final equation with numbers as requested
    print(f"ΔE = ({hbar:.5e} J·s)² / (({m_particle:.3e} kg) * ({R:.1e} m)²) / {eV_to_J:.5e} J/eV")

    # Print the final result
    print(f"\nΔE = {delta_E_eV:.5f} eV")

    return delta_E_eV

# Execute the function and save the final answer
final_answer = calculate_energy_difference()
# The final answer will be printed by the function above.
# The numeric value is captured here for the final tag.
# For example: <<<0.00847>>>