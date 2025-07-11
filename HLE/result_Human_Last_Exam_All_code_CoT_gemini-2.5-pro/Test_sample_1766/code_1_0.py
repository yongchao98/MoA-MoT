import numpy as np
from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a p-electron with spin-orbit coupling.
    """
    # Step 1: Define quantum numbers for a p-electron.
    l = 1  # Orbital angular momentum for p-shell
    s = 0.5  # Spin angular momentum for an electron

    # Step 2: Calculate the possible total angular momentum quantum numbers (j).
    # j ranges from |l-s| to l+s in integer steps.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"For a p-electron, orbital angular momentum l = {l} and spin s = {Fraction(s)}.")
    j_values_str = ', '.join([str(Fraction(j)) for j in j_values])
    print(f"The possible total angular momentum quantum numbers (j) are: {j_values_str}.\n")
    print("The common eigenvalues of J^2 and J_z are found for each |j, m_j> state.")
    print("-" * 65)

    # Step 3 & 4: Iterate through each j to find eigenvalues of J^2 and J_z.
    for j in j_values:
        j_frac = Fraction(j)
        
        # Calculate the eigenvalue of J^2: j*(j+1)ħ^2
        j_plus_1_frac = Fraction(j + 1)
        j2_eigenvalue_frac = Fraction(j * (j + 1))
        
        print(f"For the multiplet with j = {j_frac}:")
        
        # Output the J^2 eigenvalue calculation step-by-step
        print(f"  The eigenvalue of J^2 is j(j+1)ħ^2")
        print(f"  = {j_frac} * ({j_frac} + 1)ħ^2")
        print(f"  = {j_frac} * {j_plus_1_frac}ħ^2")
        print(f"  = {j2_eigenvalue_frac}ħ^2\n")
        
        # Calculate the possible m_j values and J_z eigenvalues
        mj_values = np.arange(-j, j + 1, 1)
        mj_values_str = ', '.join([str(Fraction(m)) for m in mj_values])
        jz_eigenvalues_str = ', '.join([f'{Fraction(m)}ħ' for m in mj_values])
        
        print(f"  The possible magnetic quantum numbers (m_j) are: {mj_values_str}.")
        print(f"  The corresponding eigenvalues of J_z are m_j*ħ: {jz_eigenvalues_str}.")
        print("-" * 65)

solve_spin_orbit_coupling()