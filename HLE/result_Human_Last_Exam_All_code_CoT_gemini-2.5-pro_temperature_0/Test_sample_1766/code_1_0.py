import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a p-electron (l=1, s=1/2).
    Eigenvalues are calculated in units where h_bar = 1.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1
    s = 0.5

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin s = {s}.")
    print("-" * 60)

    # Step 2: Determine possible total angular momentum quantum numbers (j)
    j_min = abs(l - s)
    j_max = l + s
    # The values of j go from j_min to j_max in steps of 1
    j_values = np.arange(j_min, j_max + 1, 1)

    print(f"The possible values for the total angular momentum quantum number j are from |l-s| to l+s.")
    print(f"j can be: {list(j_values)}")
    print("-" * 60)
    print("The common eigenvalues for J^2 and J_z are found for each state |j, m_j>.")
    print("Eigenvalue of J^2 = j * (j + 1)")
    print("Eigenvalue of J_z = m_j")
    print("-" * 60)

    # Step 3 & 4: Calculate and print eigenvalues for each j
    all_eigenvalues = []
    for j in j_values:
        # Calculate J^2 eigenvalue
        j2_eigenvalue = j * (j + 1)
        
        print(f"For the state with j = {j}:")
        print(f"  The eigenvalue of J^2 is {j} * ({j} + 1) = {j2_eigenvalue}")

        # Determine possible m_j values for the current j
        mj_values = np.arange(-j, j + 1, 1)
        
        print(f"  The possible m_j values are: {list(mj_values)}")

        # Print the common eigenvalues for each m_j
        for mj in mj_values:
            # The J_z eigenvalue is simply mj
            jz_eigenvalue = mj
            print(f"    For state |j={j}, m_j={mj}>, the (J^2, J_z) eigenvalues are ({j2_eigenvalue}, {jz_eigenvalue})")
            all_eigenvalues.append((j2_eigenvalue, jz_eigenvalue))
        print()

if __name__ == '__main__':
    solve_spin_orbital_coupling()
<<<[(3.75, 1.5), (3.75, 0.5), (3.75, -0.5), (3.75, -1.5), (0.75, 0.5), (0.75, -0.5)]>>>