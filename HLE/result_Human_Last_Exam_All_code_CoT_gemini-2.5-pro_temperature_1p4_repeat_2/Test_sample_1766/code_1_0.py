import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues for J^2 and J_z
    for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1
    s = 0.5
    
    print(f"Solving for the common eigenvalues of J^2 and J_z for a p-electron (l={l}, s={s}).")
    print("The eigenvalues of J^2 are of the form j(j+1)ħ².")
    print("The eigenvalues of J_z are of the form m_j*ħ.")
    print("-" * 60)

    # Step 2: Calculate the possible values for the total angular momentum quantum number j
    j_min = abs(l - s)
    j_max = l + s
    # j takes values from |l-s| to l+s in integer steps
    j_values = np.arange(j_min, j_max + 1, 1.0)

    # Iterate through each possible j value
    for j in j_values:
        # Step 4: Calculate J^2 eigenvalue
        j2_eigenvalue_val = j * (j + 1)
        
        print(f"For the state with total angular momentum quantum number j = {j}:")
        print(f"  The eigenvalue of J^2 is j(j+1)ħ² = {j}({j}+1)ħ² = {j2_eigenvalue_val:.2f}ħ².")
        
        # Step 3: Determine the possible m_j values for the current j
        mj_values = np.arange(-j, j + 1, 1.0)
        
        print(f"  The possible magnetic quantum numbers m_j are: {', '.join(map(str, mj_values))}.")
        print("  The corresponding common eigenvalue pairs (for J^2, J_z) are:")

        # Iterate through m_j values and print the pair of eigenvalues
        for mj in mj_values:
            # J_z eigenvalue is mj*ħ
            jz_eigenvalue_val = mj
            print(f"    - ( {j2_eigenvalue_val:.2f}ħ², {jz_eigenvalue_val:+.1f}ħ ) for the state |j, m_j> = |{j}, {mj:+.1f}>")
        print("-" * 60)

if __name__ == '__main__':
    solve_spin_orbital_coupling()