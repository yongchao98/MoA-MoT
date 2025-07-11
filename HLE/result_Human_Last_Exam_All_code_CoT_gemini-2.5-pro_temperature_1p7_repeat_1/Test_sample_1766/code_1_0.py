import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1
    s = 0.5

    print(f"For a p-electron, orbital angular momentum l = {l} and spin s = {s}.")
    print("The total angular momentum J = L + S.")
    print("The common eigenstates of J^2 and J_z are denoted by |j, m_j>.")
    print("-" * 60)

    # Step 2: Calculate the possible values for the total angular momentum quantum number j
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    # Step 3, 4, 5: Iterate through j and m_j, then calculate and print eigenvalues
    for j in j_values:
        # Calculate J^2 eigenvalue factor
        j2_eigenvalue_factor = j * (j + 1)
        
        # Use Fraction for clean output of half-integer numbers
        j_frac = Fraction(j).limit_denominator()
        j2_eigenvalue_factor_frac = Fraction(j2_eigenvalue_factor).limit_denominator()

        print(f"\nFor the states with total angular momentum j = {j_frac}:")
        
        # Determine possible m_j values
        mj_values = np.arange(-j, j + 1, 1)
        
        for mj in mj_values:
            # J_z eigenvalue factor is just mj
            mj_frac = Fraction(mj).limit_denominator()
            
            # Print the eigenvalue equations for the state |j, mj>
            print(f"\n  For the state |j={j_frac}, m_j={mj_frac}>:")
            print(f"    The J^2 eigenvalue equation is: J^2 |{j_frac}, {mj_frac}> = {j2_eigenvalue_factor_frac} \u0127\u00b2 |{j_frac}, {mj_frac}>")
            print(f"    The J_z eigenvalue equation is: J_z |{j_frac}, {mj_frac}> = {mj_frac} \u0127 |{j_frac}, {mj_frac}>")
        print("-" * 60)

if __name__ == '__main__':
    solve_spin_orbital_coupling()