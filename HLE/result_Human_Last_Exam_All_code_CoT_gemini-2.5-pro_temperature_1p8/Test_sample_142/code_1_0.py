import math

# Explanation:
# This script calculates the number of Configuration State Functions (CSFs) for a
# Full Configuration Interaction (Full CI) calculation.

# We will first attempt to use the PySCF library to automatically determine the
# number of electrons and basis functions. If PySCF is not installed, the script
# will fall back to using pre-calculated values for this specific problem.
try:
    import pyscf
except ImportError:
    print("PySCF library not found. To run this script dynamically, please install it with: pip install pyscf")
    print("Using pre-calculated values for CH2SiHF with the 6-31G** basis set.\n")
    # Manually determined values for CH2SiHF with 6-31G** basis.
    num_orbitals = 64
    num_electrons = 32
else:
    print("Using PySCF to determine system parameters.\n")
    # Define the molecule's atoms. The geometry is not needed for this calculation.
    mol = pyscf.gto.Mole()
    mol.atom = '''
        C   0.0 0.0 0.0
        H   1.0 0.0 0.0
        H  -0.5 0.86 0.0
        Si  0.0 0.0 1.8
        H   0.0 -1.0 1.8
        F  -1.0 0.0 1.8
    '''
    mol.basis = '6-31g**'
    mol.build()

    # Get the number of electrons and basis functions (orbitals) from the molecule object.
    num_electrons = mol.nelectron
    num_orbitals = mol.nao

# --- Main Calculation ---

# Step 1: Display system parameters.
print("--- Step 1: System Parameters ---")
print(f"Molecule: CH2SiHF")
print(f"Basis Set: 6-31G**")
print(f"Total number of electrons (N_elec) = {num_electrons}")
print(f"Number of basis functions (orbitals, N_orb) = {num_orbitals}")

# For a singlet state (S=0), alpha and beta electrons are equal.
num_alpha = num_electrons // 2
print(f"Number of alpha electrons (N_alpha) = {num_alpha}")
print("-" * 40)

# Step 2: Define the formula for singlet CSFs.
print("--- Step 2: Formula for Singlet CSFs (S=0) ---")
print("N_CSF = C(N_orb, N_alpha)^2 - C(N_orb, N_alpha-1) * C(N_orb, N_alpha+1)")
print("where C(n, k) is the binomial coefficient 'n choose k'.\n")
print(f"For this system: N_CSF = C({num_orbitals}, {num_alpha})^2 - C({num_orbitals}, {num_alpha-1}) * C({num_orbitals}, {num_alpha+1})")
print("-" * 40)

# Step 3: Calculate the terms and the final result.
print("--- Step 3: Calculation ---")
k = num_alpha
n = num_orbitals

# Calculate the necessary binomial coefficients ("combinations").
c_n_k = math.comb(n, k)
c_n_k_minus_1 = math.comb(n, k - 1)
c_n_k_plus_1 = math.comb(n, k + 1)

# Calculate the two main terms in the CSF formula.
term1 = c_n_k**2
term2 = c_n_k_minus_1 * c_n_k_plus_1

# Calculate the final number of CSFs.
num_csf = term1 - term2

# Display the final equation with all numbers filled in.
print("The values for the final equation are:")
print(f"C({n}, {k}) = {c_n_k}")
print(f"C({n}, {k-1}) = {c_n_k_minus_1}")
print(f"C({n}, {k+1}) = {c_n_k_plus_1}")
print("\nFinal Equation:")
print(f"N_CSF = ({c_n_k})^2 - ({c_n_k_minus_1} * {c_n_k_plus_1})")
print(f"N_CSF = {term1} - {term2}")
print("-" * 40)

print(f"The final number of CSFs is: {num_csf}")