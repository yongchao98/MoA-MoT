def solve_p_electron_spin_orbit():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a spin-orbit coupling problem involving a p-electron.
    Eigenvalues are expressed in units of hbar.
    """
    # Step 1: Define the orbital (l) and spin (s) quantum numbers.
    l = 1.0  # For a p-electron, l=1
    s = 0.5  # For an electron, s=1/2

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin s = {s}.\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number j.
    # j ranges from |l - s| to l + s in integer steps.
    j_min = abs(l - s)
    j_max = l + s
    
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1.0
        
    print(f"The possible values for the total angular momentum quantum number j are: {j_values}\n")
    print("--------------------------------------------------\n")

    # Iterate through each possible j value.
    for j in j_values:
        print(f"Analyzing the case for j = {j}:")

        # Step 3: Calculate the eigenvalue of J^2.
        j_squared_eigenvalue = j * (j + 1)
        print(f"  The eigenvalue of J^2 is given by j(j+1)hbar^2.")
        print(f"  Calculation: {j} * ({j} + 1) = {j_squared_eigenvalue}")
        print(f"  So, the J^2 eigenvalue is {j_squared_eigenvalue}*hbar^2.\n")

        # Step 4: Determine the possible values for the magnetic quantum number m_j.
        # m_j ranges from -j to +j in integer steps.
        m_j_values = []
        current_m_j = -j
        while current_m_j <= j:
            m_j_values.append(current_m_j)
            current_m_j += 1.0
        
        print(f"  For j = {j}, the possible m_j values are: {m_j_values}\n")

        # Step 5 & 6: For each state |j, m_j>, print the common eigenvalues.
        for m_j in m_j_values:
            # The eigenvalue of J_z is m_j * hbar.
            j_z_eigenvalue = m_j
            print(f"    * For the common eigenstate |j={j}, m_j={m_j}>:")
            print(f"      - Eigenvalue of J^2 = {j_squared_eigenvalue}*hbar^2")
            print(f"      - Eigenvalue of J_z = {j_z_eigenvalue}*hbar")
        
        print("\n--------------------------------------------------\n")

# Execute the function to print the solution.
solve_p_electron_spin_orbit()