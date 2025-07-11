import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    
    A p-electron has orbital angular momentum quantum number l=1 and
    spin quantum number s=1/2.
    """
    
    # Step 1: Define the quantum numbers
    l = 1.0  # Orbital angular momentum for a p-electron
    s = 0.5  # Spin angular momentum for an electron
    
    print(f"Solving for a p-electron with l = {l} and s = {s}.")
    print("The eigenvalues of J^2 are of the form j*(j+1)ħ².")
    print("The eigenvalues of J_z are of the form m_j*ħ.")
    print("-" * 60)
    
    # Step 2: Find the possible total angular momentum quantum numbers (j)
    j_values = np.arange(abs(l - s), l + s + 1, 1)
    
    # Step 3 & 4: Loop through j and m_j to find all eigenvalue pairs
    for j in j_values:
        # Calculate the eigenvalue for J^2
        j2_eigenvalue = j * (j + 1)
        
        # m_j ranges from -j to +j in steps of 1
        mj_values = np.arange(-j, j + 1, 1)
        
        for mj in mj_values:
            # The eigenvalue for J_z is simply mj
            jz_eigenvalue = mj
            
            # Step 5: Print the results for each common eigenstate |j, m_j>
            print(f"For the state with quantum numbers j = {j} and m_j = {mj}:")
            # The request is to output each number in the final equation.
            print(f"  Eigenvalue of J^2 = {j} * ({j} + 1) ħ² = {j2_eigenvalue} ħ²")
            print(f"  Eigenvalue of J_z = {jz_eigenvalue} ħ")
            print("") # Add a blank line for readability

# Execute the function
solve_spin_orbit_coupling()