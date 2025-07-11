import numpy as np

def solve_J2_Jz_eigenvalues():
    """
    Calculates and prints the common eigenvalues for J^2 and J_z for a p-electron.
    The eigenvalues are presented in the format of the eigenvalue equation.
    We are using units where hbar = 1, but will show 'hbar' symbolically in the output.
    """
    
    # Step 1: Quantum numbers for a p-electron
    l = 1
    s = 0.5
    
    print(f"For a p-electron (l={l}, s={s}), we find the common eigenvalues of J^2 and J_z.\n")

    # Step 2: Possible values for the total angular momentum quantum number j
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    # Step 3 & 4: Iterate through all states (j, m_j) and print eigenvalues
    for j in j_values:
        # Calculate the eigenvalue for J^2
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"--- States for j = {j} ---")
        
        # Iterate through possible m_j values for the current j
        m_j_values = np.arange(-j, j + 1, 1)
        
        for m_j in m_j_values:
            print(f"For the state |j={j}, m_j={m_j}>:")
            
            # Eigenvalue equation for J^2
            print(f"  J^2 |{j}, {m_j}> = {j}*({j}+1) \u0127\u00b2 |{j}, {m_j}> = {j_squared_eigenvalue} \u0127\u00b2 |{j}, {m_j}>")
            
            # Eigenvalue equation for J_z
            print(f"  J_z |{j}, {m_j}> = {m_j} \u0127 |{j}, {m_j}>\n")

# Execute the function to find and print the solutions
solve_J2_Jz_eigenvalues()