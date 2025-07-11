import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1.0  # Orbital angular momentum for p-electron
    s = 0.5  # Spin angular momentum for an electron

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin quantum number s = {s}.")
    print("The eigenvalues of J^2 are given by j(j+1)ħ^2.")
    print("The eigenvalues of J_z are given by m_j*ħ.")
    print("-" * 60)

    # Step 2: Determine possible values of j
    j_min = abs(l - s)
    j_max = l + s
    j_values = np.arange(j_min, j_max + 1, 1.0)

    # Step 3 & 4: Loop through j and m_j to calculate and print eigenvalues
    for j in j_values:
        # Calculate J^2 eigenvalue
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"For total angular momentum quantum number j = {j}:")
        
        # Print J^2 eigenvalue with its formula
        print(f"  The eigenvalue of J^2 is {j} * ({j} + 1.0) ħ^2 = {j_squared_eigenvalue} ħ^2")
        
        # Determine m_j values and calculate J_z eigenvalues
        m_j_values = np.arange(-j, j + 1, 1.0)
        print(f"  The corresponding eigenvalues of J_z are:")
        
        for m_j in m_j_values:
            # The J_z eigenvalue is just m_j
            print(f"    For m_j = {m_j:4.1f}, the eigenvalue is {m_j:4.1f} ħ")
            
        print("-" * 60)

solve_spin_orbital_coupling()