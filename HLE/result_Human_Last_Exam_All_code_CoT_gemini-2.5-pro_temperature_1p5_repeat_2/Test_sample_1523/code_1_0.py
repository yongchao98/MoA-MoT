import numpy as np
from scipy.integrate import quad

# Step 1 & 2: Define the problem and the physical quantity.
# The inner product (ϕ, D_ϕ) for a neutral scalar field is twice its Euclidean action S[ϕ].
# Its expectation value is related to the thermodynamic energy of the system.
# We calculate the finite part of the energy density, known as the thermal energy density (ε_T).

# Step 3: Define parameters for the calculation. We'll use natural units where ħ=c=k_B=1.
# Energy and momentum are measured in Mega-electron Volts (MeV).
m = 150.0  # Mass of the scalar particle in MeV (e.g., a neutral pion is ~135 MeV)
T = 200.0  # Temperature of the system in MeV

# Inverse temperature in MeV⁻¹
beta = 1.0 / T

# Step 4: Define the function to be integrated.
# This corresponds to the integrand from the formula for thermal energy density.
# The full formula is: ε_T = (1/(2π²)) * ∫[0, ∞] dk * k² * E(k) / (exp(β*E(k)) - 1)
def thermal_integrand(k, mass, beta_val):
    """
    Integrand for the thermal energy density calculation.
    k: momentum magnitude
    mass: particle mass
    beta_val: inverse temperature (1/T)
    """
    # Relativistic energy
    energy_k = np.sqrt(k**2 + mass**2)
    
    # Bose-Einstein distribution factor in the denominator.
    # np.expm1(x) calculates exp(x) - 1, which is more accurate for small x.
    bose_einstein_factor = np.expm1(beta_val * energy_k)

    # Avoid division by zero if the factor is zero (e.g., if k, m, and T are all 0)
    if bose_einstein_factor == 0:
        return 0.0
        
    return k**2 * energy_k / bose_einstein_factor

# Perform the numerical integration from k=0 to k=infinity.
# The quad function returns the integral result and an estimated error.
integral_result, error = quad(thermal_integrand, 0, np.inf, args=(m, beta))

# The final thermal energy density is the integral result multiplied by the prefactor.
thermal_energy_density = integral_result / (2 * np.pi**2)

# Print the results in a clear format
print("Problem: For a neutral scalar field at finite temperature, what is the inner product (ϕ,D_ϕ)?")
print("\nAnswer: The inner product (ϕ,D_ϕ) is equal to 2*S[ϕ], where S[ϕ] is the Euclidean action.")
print("A related physical quantity is the system's thermal energy density (ε_T), which we calculate below.")
print("-" * 70)
print("Calculation of Thermal Energy Density ε_T")
print("-" * 70)
print("Formula: ε_T = (1 / (2*π²)) * Integral from 0 to ∞ of integrand(k) dk")
print("\nParameters used:")
print(f"  Mass (m)              = {m} MeV")
print(f"  Temperature (T)       = {T} MeV")
print(f"  Inverse Temp. (β=1/T) = {beta:.4f} MeV⁻¹")
print("\nEquation with numbers plugged in:")
print(f"  integrand(k) = (k² * sqrt(k² + {m:.1f}²)) / (exp({beta:.4f} * sqrt(k² + {m:.1f}²)) - 1)")
print("-" * 70)
print("Final Result:")
# The units of energy density are Energy/Volume -> MeV / MeV⁻³ -> MeV⁴
print(f"The thermal energy density is: ε_T = {thermal_energy_density:.4e} MeV⁴")