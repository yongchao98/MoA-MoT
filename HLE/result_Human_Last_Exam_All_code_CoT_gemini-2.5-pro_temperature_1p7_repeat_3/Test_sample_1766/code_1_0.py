import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define initial quantum numbers for a p-electron.
    l = 1  # Orbital angular momentum quantum number for a p-orbital
    s = 0.5  # Spin quantum number for an electron

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin s = {s}.\n")

    # Step 2: Determine possible values for the total angular momentum quantum number, j.
    # j ranges from |l-s| to l+s in integer steps.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"The possible values for the total angular momentum quantum number, j, are: {list(j_values)}")
    print("--------------------------------------------------")

    # Step 3 & 4: For each j, find m_j values and the corresponding eigenvalues.
    for j in j_values:
        print(f"For the state with j = {j}:")

        # Calculate the eigenvalue of J^2 = j(j+1)ħ²
        j2_eigenvalue = j * (j + 1)
        print(f"  - The eigenvalue of J^2 is j(j+1)ħ²:")
        # As requested, output each number in the final equation.
        print(f"    Calculation: {j} * ({j} + 1) * ħ² = {j2_eigenvalue} ħ²")

        # Determine possible m_j values, which range from -j to +j in steps of 1
        m_j_values = np.arange(-j, j + 1, 1)

        # The eigenvalues of Jz are m_jħ
        m_j_strings = [f"{val} ħ" for val in m_j_values]
        print(f"  - The possible values for the quantum number m_j are: {list(m_j_values)}")
        print(f"  - The corresponding eigenvalues of Jz are: {', '.join(m_j_strings)}\n")
        
        print(f"  The common eigenvalue pairs (J^2 eigenvalue, Jz eigenvalue) for j={j} are:")
        for m_j in m_j_values:
            print(f"    ({j2_eigenvalue} ħ², {m_j} ħ)")
        
        print("--------------------------------------------------")

if __name__ == "__main__":
    solve_spin_orbit_coupling()