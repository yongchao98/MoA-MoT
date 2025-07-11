import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1  # Orbital angular momentum quantum number for a p-electron
    s = 0.5  # Spin quantum number for an electron

    print(f"For a p-electron, the orbital quantum number l = {l} and the spin quantum number s = {s}.\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number j
    j_min = abs(l - s)
    j_max = l + s
    # j values range from |l-s| to l+s in integer steps
    j_values = np.arange(j_min, j_max + 1, 1.0)

    print(f"The possible values for the total angular momentum quantum number j are: {list(j_values)}\n")
    print("--- Calculating Eigenvalues ---")

    # Step 3 & 4: Calculate and print the eigenvalues for each j
    for j in j_values:
        print(f"\nFor the case j = {j}:")

        # Eigenvalue of J^2
        j2_eigenvalue = j * (j + 1)
        print(f"The eigenvalue of J^2 is ħ² * j * (j + 1) = ħ² * {j} * ({j} + 1) = {j2_eigenvalue} ħ²")

        # Eigenvalues of J_z
        mj_values = np.arange(-j, j + 1, 1.0)
        print(f"The possible m_j values are {list(mj_values)}.")
        print("The corresponding eigenvalues of J_z are:")
        for mj in mj_values:
            # The eigenvalue of J_z is ħ * m_j
            print(f"  For m_j = {mj: >4.1f}, the eigenvalue is {mj: >4.1f} ħ")

solve_spin_orbital_coupling()