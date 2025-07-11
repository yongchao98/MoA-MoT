from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define the quantum numbers for a p-electron
    l = 1  # orbital angular momentum quantum number for a p-orbital
    s = Fraction(1, 2)  # spin angular momentum quantum number for an electron

    print(f"For a p-electron, the quantum numbers are l = {l} and s = {s}.")
    print("The total angular momentum J is the sum of orbital L and spin S.")
    print("\nWe need to find the eigenvalues of J^2 and Jz, which are given by j(j+1)ħ² and m_jħ respectively.")
    print("-" * 70)

    # Step 2: Calculate the possible values for the total angular momentum quantum number, j.
    j_min = abs(l - s)
    j_max = l + s
    
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1
        
    print(f"The possible values for the total angular momentum quantum number j are: {', '.join(map(str, j_values))}")

    # Step 3 & 4: For each j, find m_j values and calculate the eigenvalues.
    for j in j_values:
        j2_eigenvalue = j * (j + 1)
        
        print(f"\n--- For the state with j = {j} ---")
        
        # Print the equation and result for the J^2 eigenvalue
        print(f"The eigenvalue of J^2 is j(j+1)ħ².")
        # Output each number in the final equation
        print(f"Equation: J^2 |{j}, m_j> = {j}({j} + 1)ħ² |{j}, m_j> = {j2_eigenvalue}ħ² |{j}, m_j>")
        
        # Determine and iterate through m_j values for the current j
        mj_values = []
        current_mj = -j
        while current_mj <= j:
            mj_values.append(current_mj)
            current_mj += 1

        print(f"\nFor j = {j}, the possible m_j values are: {', '.join(map(str, mj_values))}")
        print("The common eigenvalues for each state |j, m_j> are:")

        for mj in mj_values:
            # The eigenvalue of Jz is mj*hbar
            jz_eigenvalue = mj
            print(
                f"  State |j={j}, m_j={mj}>:  Eigenvalue of J^2 = {j2_eigenvalue}ħ², "
                f"Eigenvalue of Jz = {jz_eigenvalue}ħ"
            )
            
    print("-" * 70)

if __name__ == '__main__':
    solve_spin_orbit_coupling()