import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1  # orbital angular momentum for a p-orbital
    s = 0.5  # spin angular momentum for an electron

    # Use unicode characters for clarity in the output
    hbar = "\u0127"
    hbar_sq = "\u0127\u00b2"

    print("Solving for the common eigenvalues of J\u00b2 and Jz for a p-electron.")
    print(f"The electron has orbital quantum number l = {l} and spin quantum number s = {s}.")
    print("-" * 60)

    # Step 2: Determine possible values for j
    j_min = abs(l - s)
    j_max = l + s
    # Values of j go from |l-s| to l+s in integer steps
    j_values = np.arange(j_min, j_max + 1, 1)

    print(f"The possible values for the total angular momentum quantum number j are: {list(j_values)}\n")

    # Step 3 & 4: For each j, find mj values and calculate all eigenvalues
    for j in j_values:
        # Calculate J^2 eigenvalue
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"For the state with j = {j}:")
        print(f"  The eigenvalue of J\u00b2 is {hbar_sq} * j * (j + 1)")
        print(f"  = {hbar_sq} * {j} * ({j} + 1) = {j_squared_eigenvalue} {hbar_sq}")
        
        # Determine possible m_j values
        mj_values = np.arange(-j, j + 1, 1)
        
        print(f"  The possible m_j values are {list(mj_values)}.")
        print(f"  The corresponding eigenvalues of Jz are:")
        
        # Calculate J_z eigenvalues for each m_j
        for mj in mj_values:
            jz_eigenvalue = mj
            # The format specifier ensures nice alignment for positive and negative numbers
            print(f"    For m_j = {mj: >4.1f}, the eigenvalue is {jz_eigenvalue: >4.1f} {hbar}")
            
        print("-" * 60)

solve_spin_orbital_coupling()