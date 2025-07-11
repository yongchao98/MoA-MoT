import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    
    A p-electron has orbital angular momentum quantum number l=1 and
    spin quantum number s=1/2.
    """
    l = 1
    s = 0.5

    print(f"For a p-electron, l = {l} and s = {s}.")
    print("We will find the common eigenvalues of J^2 and J_z, where J = L + S.")
    print("The eigenvalues are calculated in atomic units (h_bar = 1).\n")

    # Step 1: Find the possible values for the total angular momentum quantum number, j.
    # j ranges from |l-s| to l+s in integer steps.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print("-" * 50)
    # Step 2: For each j, find the possible m_j values and calculate eigenvalues.
    for j in j_values:
        # m_j ranges from -j to +j in integer steps.
        m_j_values = np.arange(-j, j + 1, 1)
        
        # Eigenvalue of J^2 is j(j+1)
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"For the states with total angular momentum quantum number j = {j}:")
        print(f"The eigenvalue of J^2 is j*(j+1) = {j} * ({j} + 1) = {j_squared_eigenvalue}")
        print("The corresponding eigenvalues of J_z for each m_j are:")
        
        for m_j in m_j_values:
            # Eigenvalue of J_z is m_j
            j_z_eigenvalue = m_j
            print(f"  For m_j = {m_j:<4.1f}, the J_z eigenvalue is {j_z_eigenvalue:.1f}")
        print("-" * 50)

solve_spin_orbit_coupling()