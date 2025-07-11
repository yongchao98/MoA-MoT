import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define the quantum numbers for a p-electron.
    l = 1  # orbital angular momentum quantum number for a p-orbital
    s = 0.5  # spin quantum number for an electron

    print("Problem: Find the common eigenvalues of J^2 and J_z for a p-electron.")
    print(f"A p-electron has orbital quantum number l = {l} and spin s = {s}.")
    print("The total angular momentum is J = L + S.")
    print("We are looking for states |j, m_j> that are eigenstates of both J^2 and J_z.")
    print("The eigenvalue equations are:")
    print("  J^2 |j, m_j> = j(j+1) hbar^2 |j, m_j>")
    print("  J_z |j, m_j> = m_j hbar |j, m_j>\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number, j.
    # The rule is j = |l-s|, |l-s|+1, ..., l+s
    j_values = np.arange(abs(l - s), l + s + 1, 1.0)

    # Step 3: Loop through each possible j, calculate eigenvalues, and print the results.
    for j in j_values:
        # Use fractions for clear output
        j_frac = Fraction(j).limit_denominator()
        
        print("-" * 50)
        print(f"For the case with total angular momentum j = {j_frac}:")
        
        # Calculate and print the J^2 eigenvalue
        j_squared_eigenvalue = j * (j + 1)
        j_sq_eig_frac = Fraction(j_squared_eigenvalue).limit_denominator()
        
        print("\nEigenvalue for J^2:")
        print(f"  The calculation is j(j+1) = {j_frac}({j_frac} + 1) = {j_sq_eig_frac}.")
        print(f"  The eigenvalue equation is: J^2 |j={j_frac}, m_j> = {j_sq_eig_frac} hbar^2 |j={j_frac}, m_j>")
        
        # Determine the possible m_j values for the current j
        m_j_values = np.arange(-j, j + 1, 1.0)
        
        print("\nEigenvalues for J_z:")
        # Loop through m_j values and print the eigenvalue equation for each
        for m_j in m_j_values:
            m_j_frac = Fraction(m_j).limit_denominator()
            print(f"  For m_j = {m_j_frac}, the eigenvalue equation is: J_z |j={j_frac}, m_j={m_j_frac}> = {m_j_frac} hbar |j={j_frac}, m_j={m_j_frac}>")

solve_spin_orbital_coupling()