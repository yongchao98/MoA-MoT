from fractions import Fraction

def solve_p_electron_eigenvalues():
    """
    Solves for the common eigenvalues of J^2 and J_z for a p-electron.

    A p-electron has orbital angular momentum quantum number l=1 and
    spin quantum number s=1/2. The total angular momentum J is the
    sum of the orbital (L) and spin (S) angular momenta.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1
    s = Fraction(1, 2)

    print(f"Solving for a p-electron with l = {l} and s = {s}.")
    print("The eigenvalues are given in units of hbar (ħ).")
    print("-" * 50)

    # Step 2: Determine possible values for j
    j_min = abs(l - s)
    j_max = l + s
    # j ranges from |l-s| to l+s in steps of 1
    j_values = [j_min + i for i in range(int(j_max - j_min) + 1)]

    # Iterate through each possible total angular momentum state
    for j in j_values:
        print(f"For the state with total angular momentum quantum number j = {j}:")

        # Step 3: Calculate and print the eigenvalue for J^2
        j_squared_eigenvalue = j * (j + 1)
        print(f"  - The eigenvalue of J^2 is j(j+1)ħ²")
        print(f"    j(j+1)ħ² = {j} * ({j} + 1)ħ² = {j_squared_eigenvalue}ħ²")

        # Step 4: Determine the possible values for m_j
        m_j_values = []
        current_m_j = -j
        while current_m_j <= j:
            m_j_values.append(current_m_j)
            current_m_j += 1

        # Print the possible eigenvalues for J_z
        m_j_str_list = [str(val) for val in m_j_values]
        print(f"  - The eigenvalues of J_z are m_jħ, where m_j can be:")
        print(f"    m_j = {', '.join(m_j_str_list)}")
        print("-" * 50)

if __name__ == '__main__':
    solve_p_electron_eigenvalues()