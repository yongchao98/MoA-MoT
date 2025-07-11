from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and Jz
    for a p-electron (l=1, s=1/2).
    """
    # Step 1: Define the quantum numbers for a p-electron.
    l = 1
    s = Fraction(1, 2)

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin quantum number s = {s}.")
    
    # Step 2: Determine the possible values for the total angular momentum quantum number j.
    # j ranges from |l-s| to l+s in integer steps.
    j_min = abs(l - s)
    j_max = l + s
    
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1
        
    print(f"The possible values for the total angular momentum quantum number j are: {', '.join(map(str, j_values))}\n")
    print("-----------------------------------------------------------")
    print("The common eigenvalues of J^2 and J_z are:")
    print("-----------------------------------------------------------")

    # Iterate through each possible j value
    for j in j_values:
        # Step 3: Calculate the eigenvalue of J^2.
        # Eigenvalue of J^2 is j(j+1)ħ^2.
        j_squared_eigenvalue = j * (j + 1)
        
        # Step 4: Determine the possible values for the magnetic quantum number m_j.
        # m_j ranges from -j to +j in integer steps.
        m_j = -j
        while m_j <= j:
            # Step 5: The eigenvalue of J_z is m_j*ħ.
            jz_eigenvalue = m_j
            
            # Print the results in a detailed format
            print(f"For quantum numbers (j = {j}, m_j = {m_j}):")
            
            # Show the calculation for the J^2 eigenvalue
            print(f"  Eigenvalue of J^2 = j(j+1)ħ² = {j} * ({j} + 1)ħ² = {j_squared_eigenvalue}ħ²")
            
            # Show the calculation for the J_z eigenvalue
            print(f"  Eigenvalue of J_z = m_j*ħ = {jz_eigenvalue}ħ\n")
            
            m_j += 1

if __name__ == "__main__":
    solve_spin_orbit_coupling()
