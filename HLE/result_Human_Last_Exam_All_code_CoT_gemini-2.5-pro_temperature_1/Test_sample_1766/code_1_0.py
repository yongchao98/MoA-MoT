import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1  # Orbital angular momentum for p-shell
    s = 0.5  # Spin angular momentum for an electron

    print(f"Solving for a p-electron with l = {l} and s = {s}\n")

    # Step 2: Determine possible total angular momentum quantum numbers (j)
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"The possible values for the total angular momentum quantum number j are {list(j_values)}.\n")

    # Header for the results
    print("The common eigenvalues for J^2 and J_z are determined for each state |j, mj>:")
    print("-" * 70)

    # Sort j values in descending order for conventional representation
    j_values = sorted(j_values, reverse=True)
    
    # Step 3 & 4: Calculate and print eigenvalues for each j and corresponding m_j
    for j in j_values:
        # Calculate J^2 eigenvalue
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"For the state with j = {j}:")
        print(f"  The J^2 eigenvalue is j(j+1)hbar^2 = {j} * ({j} + 1) hbar^2 = {j_squared_eigenvalue:.2f} hbar^2")

        # Determine possible m_j values
        m_j_values = np.arange(-j, j + 1, 1)
        
        print(f"  The corresponding m_j values range from -{j} to +{j}:")
        
        # Calculate J_z eigenvalue for each m_j
        for m_j in m_j_values:
            j_z_eigenvalue = m_j
            print(f"    For m_j = {m_j: >4.1f}, the J_z eigenvalue is {j_z_eigenvalue: >4.1f} hbar")
        
        print("-" * 70)

if __name__ == "__main__":
    solve_spin_orbit_coupling()