from fractions import Fraction

def solve_p_electron_eigenvalues():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    The final equations for each eigenstate |j, m_j> are displayed.
    """
    # For a p-electron, orbital angular momentum quantum number l = 1.
    l = 1
    # For an electron, spin quantum number s = 1/2.
    s = 0.5

    # Determine the possible values for the total angular momentum quantum number, j.
    # j ranges from |l-s| to l+s in integer steps.
    j_values = []
    j = abs(l - s)
    while j <= l + s:
        j_values.append(j)
        j += 1

    print(f"For a p-electron with l={l} and s={Fraction(s)}, the possible values for the total angular momentum quantum number j are: {[str(Fraction(val)) for val in j_values]}\n")

    hbar_str = "ħ"
    hbar2_str = "ħ²"

    # Iterate through each possible j value
    for j in j_values:
        j_frac = Fraction(j)
        # Calculate the eigenvalue for J^2, which is j(j+1)
        j_squared_val = j * (j + 1)
        j_squared_frac = Fraction(j_squared_val)

        print(f"--- For states with total angular momentum j = {j_frac} ---")
        
        # Iterate through possible m_j values for the current j
        m_j = -j
        while m_j <= j:
            m_j_frac = Fraction(m_j)
            
            # Print the state being considered
            state_ket = f"|j={j_frac}, m_j={m_j_frac}⟩"
            print(f"For the common eigenstate {state_ket}:")
            
            # Print J^2 eigenvalue equation
            print(f"  J² {state_ket} = {j_frac}*({j_frac} + 1) {hbar2_str} {state_ket} = {j_squared_frac} {hbar2_str} {state_ket}")

            # Print J_z eigenvalue equation
            print(f"  Jz {state_ket} = {m_j_frac} {hbar_str} {state_ket}\n")

            m_j += 1

# Execute the function to print the solution
solve_p_electron_eigenvalues()