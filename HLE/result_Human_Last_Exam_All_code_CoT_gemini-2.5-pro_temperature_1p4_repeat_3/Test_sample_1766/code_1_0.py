import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define the quantum numbers for a p-electron
    l = 1  # orbital angular momentum quantum number
    s = 0.5  # spin angular momentum quantum number

    print(f"For a p-electron, the quantum numbers are l = {l} and s = {s}.")
    print("-" * 55)

    # Step 2: Determine the possible total angular momentum quantum numbers, j
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print("The possible values for the total angular momentum quantum number j are:")
    for j_val in j_values:
        print(f"j = {j_val}")
    print("-" * 55)

    print("The common eigenvalues of J^2 and J_z are found for each j-multiplet:")
    print("-" * 55)

    # Step 3 & 4: Calculate and print the eigenvalues for each multiplet
    for j in j_values:
        # Calculate the J^2 eigenvalue
        j_squared_eigenvalue_num = j * (j + 1)
        
        print(f"For the j = {j} multiplet:")
        
        # Format the J^2 eigenvalue as a fraction for clean output
        j2_frac = Fraction(j_squared_eigenvalue_num).limit_denominator()
        j2_str = f"({j2_frac.numerator}/{j2_frac.denominator})ħ²"
        
        # Output the J^2 eigenvalue including the calculation
        print(f"  The eigenvalue of J^2 is j(j+1)ħ² = {j}({j}+1.0)ħ² = {j2_str}")
        
        # Calculate the corresponding J_z eigenvalues
        mj_values = np.arange(-j, j + 1, 1)
        
        print(f"  The corresponding eigenvalues of Jz (m_j * ħ) are:")
        for mj in mj_values:
            # Format the J_z eigenvalue as a fraction
            mj_frac = Fraction(mj).limit_denominator()
            jz_str = f"({mj_frac.numerator}/{mj_frac.denominator})ħ"
            
            # Output the J_z eigenvalue for each m_j
            print(f"    For m_j = {mj: >4.1f}, the eigenvalue is {jz_str}")
        print("-" * 55)

# Execute the function
solve_spin_orbital_coupling()