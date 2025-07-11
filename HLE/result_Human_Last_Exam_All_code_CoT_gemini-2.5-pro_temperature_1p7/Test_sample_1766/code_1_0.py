import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a p-electron.
    """
    # For a p-electron, the orbital angular momentum quantum number l is 1.
    l = 1
    # For an electron, the spin quantum number s is 1/2.
    s = 0.5

    print("For a p-electron, we have orbital quantum number l = 1 and spin quantum number s = 1/2.")
    print("The total angular momentum J is the sum of orbital (L) and spin (S) angular momentum.")
    print("The possible values for the total angular momentum quantum number, j, range from |l-s| to l+s.")
    print(f"So, for l=1 and s=0.5, j can be |{l} - {s}| = {abs(l-s)} or {l} + {s} = {l+s}.\n")

    # The possible values for the total angular momentum quantum number j
    j_values = [l + s, abs(l - s)]

    print("The eigenvalues of the J^2 operator are j(j+1)*hbar^2.")
    print("The eigenvalues of the J_z operator are m_j*hbar.")
    print("We will calculate these eigenvalues for each possible j (in units where hbar = 1).\n")
    print("-" * 70)

    # Loop through each possible j value, starting from the largest
    for j in sorted(j_values, reverse=True):
        # Calculate the eigenvalue for J^2
        j2_eigenvalue = j * (j + 1)
        
        # Determine the possible m_j values for the current j
        # m_j ranges from -j to +j in integer steps.
        m_j_values = np.arange(-j, j + 1, 1.0)
        
        print(f"For the state with total angular momentum quantum number j = {j}:")
        
        # Print the J^2 eigenvalue equation and result
        print(f"  The eigenvalue of J^2 = j(j+1) = {j} * ({j} + 1) = {j2_eigenvalue:.2f}")
        
        # Print the J_z eigenvalues for all possible m_j
        print(f"  For this j, the magnetic quantum number m_j can be: {', '.join(map(str, m_j_values))}.")
        print("  The corresponding eigenvalues of J_z are m_j:")
        for m_j in m_j_values:
            print(f"    m_j = {m_j: >4.1f}  =>  Eigenvalue of J_z = {m_j: >4.1f}")
            
        print("-" * 70)

# Run the function to display the solution
solve_spin_orbit_coupling()