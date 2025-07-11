from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # For a p-electron, the orbital angular momentum quantum number l=1.
    l = 1
    # For an electron, the spin quantum number s=1/2.
    s = Fraction(1, 2)

    print(f"Solving for a p-electron with l={l} and s={s}:\n")

    # The total angular momentum quantum number j ranges from |l-s| to l+s.
    j_min = abs(l - s)
    j_max = l + s
    
    # Generate the possible j values (with integer steps)
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1
    
    # For each state |j, m_j>, calculate the eigenvalues of J^2 and J_z.
    for j in j_values:
        j2_eigenvalue = j * (j + 1)
        
        print(f"--- For total angular momentum j = {j} ---")
        print(f"The J^2 eigenvalue is j(j+1)ħ^2 = {j}({j}+1)ħ^2 = {j2_eigenvalue}ħ^2.")
        print(f"The common eigenstates |j, m_j> and J_z eigenvalues are:")

        # The quantum number m_j ranges from -j to +j in integer steps.
        current_mj = -j
        while current_mj <= j:
            # The J_z eigenvalue is m_j * hbar
            print(f"  For state |j={j}, m_j={current_mj}>, the J_z eigenvalue is {current_mj}ħ.")
            current_mj += 1
        print("") # Add a blank line for readability

solve_spin_orbit_coupling()
<<<The possible total angular momentum quantum numbers are j=3/2 and j=1/2. For j=3/2, m_j can be 3/2, 1/2, -1/2, -3/2. For j=1/2, m_j can be 1/2, -1/2.>>>