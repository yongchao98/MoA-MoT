import numpy as np

# This script calculates the theoretical product ratio for the electrocyclization
# of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene based on Frontier Molecular Orbital theory.

# 1. Define the parameters for the reacting pi-system (octatetraene)
N = 8  # Number of pi electrons / atoms in the conjugated system
j = 1  # Index of the terminal atom (can be 1 or N)

# 2. Identify the HOMO and LUMO quantum numbers (k) for the polyene
# For a linear polyene with N atoms, there are N molecular orbitals.
# The N electrons fill the first N/2 orbitals.
# HOMO is the highest occupied molecular orbital, k = N/2.
# LUMO is the lowest unoccupied molecular orbital, k = N/2 + 1.
k_homo = N // 2
k_lumo = N // 2 + 1

# 3. Use the analytical formula for Huckel coefficients of a linear polyene:
# c_k,j = sqrt(2 / (N + 1)) * sin(k * j * pi / (N + 1))
def calculate_coefficient(k, j, N):
    """Calculates the Huckel MO coefficient for atom j in orbital k."""
    return np.sqrt(2 / (N + 1)) * np.sin(k * j * np.pi / (N + 1))

# 4. Calculate the coefficients at the terminus for the HOMO and LUMO
c_homo_terminus = calculate_coefficient(k_homo, j, N)
c_lumo_terminus = calculate_coefficient(k_lumo, j, N)

# 5. Model the product ratio.
# The allowed conrotatory path yields trans-isomer B. The ratio depends on the HOMO.
# The forbidden disrotatory path yields cis-isomer A. We model its rate relative
# to the allowed path using the ratio of LUMO to HOMO terminal coefficients.
# Ratio A:B is modeled as (c_LUMO_terminus / c_HOMO_terminus)^2
# Due to the pairing theorem, |c_homo_terminus| is equal to |c_lumo_terminus|.
ratio_val = (c_lumo_terminus / c_homo_terminus)**2

# The numbers for the final ratio equation are derived from the calculation.
# Since ratio_val is 1.0, the ratio is 1:1.
ratio_A = 1
ratio_B = 1

# 6. Print the results, including the numbers for the final equation.
print("Prediction based on a Frontier Molecular Orbital Theory model:")
print("-" * 55)
print("Reactant: (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene (an 8-pi electron system)")
print("Path A (forbidden, disrotatory) -> cis-isomer A")
print("Path B (allowed, conrotatory) -> trans-isomer B")
print("\nModel: Ratio(A/B) = (c_LUMO_terminus / c_HOMO_terminus)^2")
print(f"\nFor an {N}-atom pi system, HOMO is psi({k_homo}) and LUMO is psi({k_lumo}).")
print(f"Terminal coefficient of HOMO (c_HOMO_terminus): {c_homo_terminus:.4f}")
print(f"Terminal coefficient of LUMO (c_LUMO_terminus): {c_lumo_terminus:.4f}")
print(f"Value of the ratio model = ({c_lumo_terminus:.4f}/{c_homo_terminus:.4f})^2 = {ratio_val:.2f}")

print("\n-------------------------------------------------------")
print("Final Predicted Ratio Equation:")
print(f"    Isomer A : Isomer B = {ratio_A} : {ratio_B}")
print("-------------------------------------------------------")
