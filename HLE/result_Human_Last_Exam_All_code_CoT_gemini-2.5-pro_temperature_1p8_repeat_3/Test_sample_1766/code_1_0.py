import fractions

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1
    s = fractions.Fraction(1, 2)

    # Step 2: Determine the possible values for the total angular momentum quantum number j
    j_min = abs(l - s)
    j_max = l + s
    # For l=1 and s=1/2, j can be 1/2 and 3/2
    j_values = [j for j in (j_min, j_max)]

    print("Solving for the common eigenvalues of J^2 and Jz for a p-electron (l=1, s=1/2).\n")
    print("The eigenvalues of J^2 are of the form j(j+1)hbar^2.")
    print("The eigenvalues of Jz are of the form m_j*hbar, where m_j = -j, -j+1, ..., +j.\n")
    print("-" * 50 + "\n")

    # Step 3: Calculate and print the eigenvalues for each possible j
    all_j2_eigenvalues = []
    all_mj_eigenvalues = set()

    for j in j_values:
        # Calculate J^2 eigenvalue
        j_plus_1 = j + 1
        j2_eigenvalue = j * j_plus_1
        all_j2_eigenvalues.append(j2_eigenvalue)

        # Determine possible m_j values
        mj_values = []
        current_mj = -j
        while current_mj <= j:
            mj_values.append(current_mj)
            all_mj_eigenvalues.add(current_mj)
            current_mj += 1

        mj_values_str = ", ".join(map(str, mj_values))

        # Print the results for the current j value
        print(f"For total angular momentum quantum number j = {j}:")
        
        # Print the J^2 eigenvalue calculation
        print(f"  The eigenvalue equation for J^2 is: j*(j+1) * hbar^2")
        print(f"  Substituting j = {j}, the eigenvalue is:")
        print(f"    {j} * ({j} + 1) * hbar^2 = {j} * {j_plus_1} * hbar^2 = {j2_eigenvalue} * hbar^2")

        # Print the J_z eigenvalues
        print(f"  The corresponding eigenvalues of Jz are m_j * hbar, where m_j can be:")
        print(f"    {{{mj_values_str}}}\n")

if __name__ == "__main__":
    solve_spin_orbit_coupling()