def solve_p_electron_eigenvalues():
    """
    Calculates and prints the common eigenvalues of J^2 and Jz
    for a single p-electron with spin-orbit coupling.
    """
    # Step 1: Define the primary quantum numbers for a p-electron.
    l = 1.0  # Orbital angular momentum quantum number for a p-orbital
    s = 0.5  # Spin quantum number for an electron

    print("Solving for the common eigenvalues of J^2 and J_z for a p-electron.")
    print(f"The quantum numbers are l = {l} (orbital) and s = {s} (spin).\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number, j.
    # The rule is j ranges from |l-s| to l+s in integer steps.
    j_min = abs(l - s)
    j_max = l + s
    
    j_values = []
    j_current = j_min
    while j_current <= j_max:
        j_values.append(j_current)
        j_current += 1.0

    print(f"The possible values for the total angular momentum quantum number j are: {j_values}\n")
    print("-------------------------------------------------------------------------")
    print("The common eigenvalues for the quantum states |j, mj> are as follows:")
    print("-------------------------------------------------------------------------")

    # Step 3 & 4: For each j, calculate J^2 and J_z eigenvalues.
    for j in j_values:
        # Calculate the eigenvalue for J^2, which is j(j+1)ħ^2
        j2_eigenvalue = j * (j + 1)
        
        # Determine the possible values for the magnetic quantum number, mj.
        # mj ranges from -j to +j in integer steps.
        mj_values = []
        mj_current = -j
        while mj_current <= j:
            mj_values.append(mj_current)
            mj_current += 1.0

        print(f"\nFor states with j = {j}:")
        print(f" The eigenvalue of J^2 is {j} * ({j} + 1.0) * ħ^2 = {j2_eigenvalue}ħ^2")
        print(f" The possible J_z eigenvalues (mj * ħ) are:")
        for mj in mj_values:
            # The J_z eigenvalue is mj*ħ
            print(f"  - For mj = {mj: >4.1f}, the eigenvalue pair is (J^2 = {j2_eigenvalue}ħ^2, J_z = {mj: >4.1f}ħ)")

# Execute the function to print the results.
solve_p_electron_eigenvalues()