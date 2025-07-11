import numpy as np
from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    The eigenvalues are expressed in terms of ħ (hbar).
    """
    l = 1  # Orbital angular momentum for a p-electron
    s = Fraction(1, 2)  # Spin angular momentum for an electron

    print(f"For a p-electron, the quantum numbers are l={l} and s={s}.")
    print("The total angular momentum J is the sum of orbital (L) and spin (S) angular momentum.")
    print("The possible values for the total angular momentum quantum number, j, range from |l-s| to l+s.")
    print("-" * 70)

    # Calculate the possible values for j
    j_min = abs(l - s)
    j_max = l + s
    j_values = [Fraction(j) for j in np.arange(j_min, j_max + 1, 1)]

    # Iterate through each possible j value
    for j in j_values:
        print(f"For the state with total angular momentum quantum number j = {j}:")

        # Calculate and display the eigenvalue of J^2
        j_plus_1 = j + 1
        j_squared_eigenvalue = j * j_plus_1
        print(f"  - The eigenvalue of J^2 is j(j+1)ħ² = {j} * ({j} + 1)ħ² = {j} * {j_plus_1}ħ² = {j_squared_eigenvalue}ħ²")

        # Calculate the possible values for m_j
        m_j_values = [Fraction(m) for m in np.arange(-j, j + 1, 1)]

        print(f"  - The corresponding eigenvalues of J_z, given by m_j*ħ, are:")
        for m_j in m_j_values:
            # Display the eigenvalue of J_z
            print(f"    For m_j = {m_j:>4}, the eigenvalue is {m_j}ħ")
        print("-" * 70)

solve_spin_orbit_coupling()
<<<For j = 3/2, the J^2 eigenvalue is 15/4*ħ² and the J_z eigenvalues are 3/2*ħ, 1/2*ħ, -1/2*ħ, -3/2*ħ. For j = 1/2, the J^2 eigenvalue is 3/4*ħ² and the J_z eigenvalues are 1/2*ħ, -1/2*ħ.>>>