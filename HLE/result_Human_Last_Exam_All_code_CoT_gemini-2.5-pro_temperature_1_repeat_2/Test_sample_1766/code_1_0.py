import numpy as np
from fractions import Fraction

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # Step 1: Define the quantum numbers for a p-electron
    l = 1  # orbital angular momentum for p-shell
    s = Fraction(1, 2)  # spin angular momentum for an electron

    print(f"For a p-electron, l = {l} and s = {s}.\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number j
    j_min = abs(l - s)
    j_max = l + s
    
    # j ranges from |l-s| to l+s in steps of 1
    j_values = np.arange(float(j_min), float(j_max) + 1, 1.0)
    
    print("The possible values for the total angular momentum quantum number j are:")
    for j_val in j_values:
        print(f"j = {Fraction(j_val)}")
    print("-" * 40)

    # Step 3 & 4: For each j, find mj values and calculate the eigenvalues
    for j_float in j_values:
        j = Fraction(j_float)
        
        # Calculate the J^2 eigenvalue
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"For the state with j = {j}:")
        print(f"  The eigenvalue of J^2 is j(j+1)hbar^2 = {j}({j}+1)hbar^2 = {j_squared_eigenvalue} hbar^2")
        
        # Determine the possible mj values
        mj_values = [Fraction(m) for m in np.arange(-float(j), float(j) + 1, 1.0)]
        
        print(f"  The possible m_j values range from -j to +j:")
        for mj in mj_values:
            # The J_z eigenvalue is simply mj
            jz_eigenvalue = mj
            print(f"    For m_j = {mj:<5}, the eigenvalue of J_z is {jz_eigenvalue} hbar")
        print("-" * 40)

solve_spin_orbital_coupling()