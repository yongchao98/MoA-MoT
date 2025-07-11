def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    A p-electron has orbital angular momentum l=1 and spin s=1/2.
    """
    l = 1
    s = 0.5

    # Determine the possible values for the total angular momentum quantum number j
    j_values = [abs(l - s), l + s]

    print(f"For a p-electron (l={l}, s={s}):")
    print(f"The possible total angular momentum quantum numbers are j = {j_values[0]} and j = {j_values[1]}.\n")
    
    # Iterate through each possible j value in descending order
    for j in sorted(j_values, reverse=True):
        print(f"--- For j = {j} ---")

        # Calculate and print the eigenvalue of J^2
        j_squared_eigenvalue = j * (j + 1)
        print(f"The eigenvalue of J^2 is j(j+1)ħ^2 = {j} * ({j} + 1)ħ^2 = {j_squared_eigenvalue}ħ^2")

        # Determine the possible m_j values for the current j
        m_j_values = []
        m_j = -j
        while m_j <= j:
            m_j_values.append(m_j)
            m_j += 1
        
        print(f"The corresponding eigenvalues of J_z for m_j from {-j} to {j} are:")
        # Iterate through possible m_j values and print J_z eigenvalues
        for mj in sorted(m_j_values, reverse=True):
            print(f"  m_j = {mj: >4.1f}, J_z eigenvalue = {mj: >4.1f}ħ")
        print() # Add a newline for better separation

if __name__ == '__main__':
    solve_spin_orbital_coupling()