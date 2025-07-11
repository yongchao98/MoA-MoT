import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    J is the total angular momentum, where J = L + S.
    """
    
    # Step 1: Define the quantum numbers for a p-electron.
    l = 1  # Orbital angular momentum quantum number for a 'p' orbital.
    s = 0.5  # Spin angular momentum quantum number for an electron.

    print(f"Solving for the common eigenvalues of J^2 and J_z for a p-electron.")
    print(f"A p-electron has orbital quantum number l = {l} and spin quantum number s = {Fraction(s)}.\n")

    # Step 2: Calculate the possible values for the total angular momentum quantum number, j.
    # The rule is |l-s| <= j <= l+s, in integer steps.
    j_min = abs(l - s)
    j_max = l + s
    j_values = np.arange(j_min, j_max + 1, 1.0) # Step is always 1
    
    j_frac_values = [Fraction(j).limit_denominator() for j in j_values]
    print(f"The possible values for the total angular momentum quantum number j are: {', '.join(map(str, j_frac_values))}\n")

    # Step 3 & 4: For each j, calculate the eigenvalues of J^2 and J_z.
    for j in j_values:
        j_frac = Fraction(j).limit_denominator()
        
        print(f"--- For the state with j = {j_frac} ---")
        
        # Calculate J^2 eigenvalue: hbar^2 * j * (j+1)
        j_squared_eigenvalue = j * (j + 1)
        j_squared_eigenvalue_frac = Fraction(j_squared_eigenvalue).limit_denominator()
        
        print(f"The eigenvalue of J^2 is calculated as j(j+1)hbar^2.")
        print(f"Result: {j_frac}({j_frac} + 1)hbar^2 = {j_frac}({j_frac + 1})hbar^2 = {j_squared_eigenvalue_frac} hbar^2")
        
        # Calculate J_z eigenvalues: hbar * m_j
        # where m_j ranges from -j to +j in integer steps.
        mj_values = np.arange(-j, j + 1, 1.0)
        mj_frac_list = [Fraction(mj).limit_denominator() for mj in mj_values]
        
        print(f"The eigenvalues of J_z are calculated as m_j*hbar, where m_j can be {', '.join(map(str, mj_frac_list))}.")
        
        # Create the string for the list of Jz eigenvalues
        jz_eigenvalues_str = ", ".join([f"{mj_f} hbar" for mj_f in mj_frac_list])
        print(f"Result: The set of J_z eigenvalues is {{{jz_eigenvalues_str}}}\n")

if __name__ == '__main__':
    solve_spin_orbital_coupling()