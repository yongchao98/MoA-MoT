import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1  # Orbital angular momentum for a p-orbital
    s = 0.5  # Spin angular momentum for an electron

    print(f"Solving for a p-electron with l={l} and s={s}\n")

    # Step 2: Determine possible values for j
    j_min = abs(l - s)
    j_max = l + s
    # j ranges from j_min to j_max in steps of 1
    j_values = np.arange(j_min, j_max + 1, 1)

    # Loop through each possible j value
    for j in j_values:
        # Step 3: For each j, determine the possible m_j values
        # m_j ranges from -j to +j in steps of 1
        mj_values = np.arange(-j, j + 1, 1)

        print(f"--- For states with total angular momentum j = {j} ---")
        
        # Step 4: Calculate the J^2 eigenvalue (it's the same for all m_j with the same j)
        j2_eigenvalue = j * (j + 1)
        
        # Loop through each m_j for the current j
        for mj in mj_values:
            # The J_z eigenvalue is just m_j
            jz_eigenvalue = mj
            
            # Print the final eigenvalue equations for the state |j, m_j>
            print(f"For the state |j={j}, m_j={mj}>:")
            print(f"  Eigenvalue of J^2 = j(j+1) = {j} * ({j} + 1) = {j2_eigenvalue}")
            print(f"  Eigenvalue of J_z = m_j = {jz_eigenvalue}\n")

solve_spin_orbit_coupling()
