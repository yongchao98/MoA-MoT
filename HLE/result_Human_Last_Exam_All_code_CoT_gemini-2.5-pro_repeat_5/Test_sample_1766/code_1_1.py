import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1      # orbital quantum number for p-shell
    s = 0.5    # spin quantum number for an electron

    print(f"Solving for the common eigenvalues of J^2 and Jz for a p-electron (l={l}, s={s}).")
    print("The eigenvalues are expressed in terms of ħ (reduced Planck constant).")
    print("-" * 70)

    # Step 2: Determine possible total angular momentum quantum numbers (j)
    j_min = abs(l - s)
    j_max = l + s
    j_values = np.arange(j_min, j_max + 1, 1)

    print(f"The possible values for the total angular momentum quantum number j are: {list(j_values)}")

    # Steps 3 & 4: Loop through each j to find eigenvalues for J^2 and J_z
    for j in j_values:
        # Calculate the eigenvalue for J^2: j*(j+1) * ħ^2
        eigenvalue_J2 = j * (j + 1)
        
        print(f"\n--- For the state with j = {j} ---")
        
        # Print the eigenvalue equation for J^2 with the numbers filled in
        print(f"The eigenvalue equation for J^2 is:")
        print(f"  J^2 |j={j}, m_j> = {j}({j}+1)ħ^2 |j={j}, m_j> = {eigenvalue_J2}ħ^2 |j={j}, m_j>")

        # Calculate the possible magnetic quantum numbers (m_j) for the given j
        mj_values = np.arange(-j, j + 1, 1)
        
        # Print the eigenvalue equations for J_z for each possible m_j
        print(f"The eigenvalue equations for J_z are:")
        for mj in mj_values:
            # The eigenvalue for J_z is mj * ħ
            eigenvalue_Jz = mj
            print(f"  J_z |j={j}, m_j={mj}> = {eigenvalue_Jz}ħ |j={j}, m_j={mj}>")

if __name__ == '__main__':
    solve_spin_orbit_coupling()