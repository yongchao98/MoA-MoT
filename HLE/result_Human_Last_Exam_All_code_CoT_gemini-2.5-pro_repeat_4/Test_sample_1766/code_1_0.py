import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues for J^2 and J_z
    for a p-electron subject to spin-orbit coupling.
    """
    # Step 1: Define the quantum numbers for a p-electron
    l = 1  # Orbital angular momentum for a p-electron
    s = 0.5  # Spin angular momentum for an electron

    print(f"Solving for the common eigenvalues of J^2 and J_z for a p-electron.")
    print(f"A p-electron has orbital angular momentum quantum number l = {l}.")
    print(f"An electron has spin angular momentum quantum number s = {Fraction(s)}.\n")

    # Step 2: Determine possible j values
    # The total angular momentum quantum number, j, can take values from |l-s| to l+s.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"The total angular momentum quantum number, j, can take values from |l-s| to l+s.")
    print(f"This results in possible j values of: {', '.join([str(Fraction(j)) for j in j_values])}\n")

    # Step 3 & 4: Calculate and print eigenvalues for each j
    for j_val in j_values:
        j = Fraction(j_val)
        
        # Calculate J^2 eigenvalue
        j_plus_1 = j + 1
        j2_eigenvalue = j * j_plus_1
        
        print(f"--- For the state with j = {j} ---")
        print(f"The eigenvalue of J^2 is calculated as j(j+1).")
        print(f"Eigenvalue of J^2 = {j} * ({j} + 1) = {j} * {j_plus_1} = {j2_eigenvalue} (in units of h_bar^2)")
        
        # Determine m_j values
        m_j_values = np.arange(-j_val, j_val + 1, 1)
        
        print(f"For j = {j}, the magnetic quantum number m_j can be: {', '.join([str(Fraction(m)) for m in m_j_values])}.")
        print(f"The eigenvalue of J_z is simply m_j (in units of h_bar).")
        
        print("\nThe common eigenvalue pairs (J^2 eigenvalue, J_z eigenvalue) are:")
        for m_j_val in m_j_values:
            m_j = Fraction(m_j_val)
            print(f"  For m_j = {str(m_j).rjust(4)}:  ({str(j2_eigenvalue).rjust(4)}, {str(m_j).rjust(4)})")
        print("-" * 40 + "\n")

solve_spin_orbital_coupling()