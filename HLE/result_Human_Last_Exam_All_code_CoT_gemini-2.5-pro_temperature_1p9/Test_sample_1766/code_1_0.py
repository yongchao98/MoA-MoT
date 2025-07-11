import numpy as np
from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    A p-electron has orbital angular momentum l=1 and spin s=1/2.
    """

    l = 1
    s = Fraction(1, 2)

    # Determine the possible values for the total angular momentum quantum number j
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"For a p-electron with l={l} and s={s}:\n")

    # Iterate through each possible j value
    for j_val in j_values:
        j = Fraction(j_val)
        print("--------------------------------------------------")
        print(f"For the state with total angular momentum j = {j}:")
        
        # Calculate and print the eigenvalue of J^2
        j2_eigenvalue = j * (j + 1)
        print(f"  The eigenvalue of J^2 is j(j+1)ħ² = {j}*({j}+1)ħ² = {j2_eigenvalue}ħ²")
        
        # Determine the possible mj values
        mj_values = [Fraction(m) for m in np.arange(-j, j + 1, 1)]
        
        # Print the corresponding eigenvalues of Jz
        print(f"  The {2*j+1} possible eigenvalues of J_z are m_j*ħ:")
        for mj in mj_values:
            print(f"    For m_j = {mj}, the eigenvalue is {mj}ħ")

solve_spin_orbit_coupling()

# The total number of states is (2l+1)*(2s+1) = (2*1+1)*(2*1/2+1) = 3 * 2 = 6.
# This corresponds to the states for j=3/2 (2*(3/2)+1 = 4 states) and j=1/2 (2*(1/2)+1 = 2 states). 4 + 2 = 6.
total_states = (2 * l + 1) * (2 * s.numerator / s.denominator + 1)

# As the problem does not have a single numerical answer, 
# we provide the total number of eigenstates as the final answer.
final_answer = int((2 * 1 + 1) * (2 * 0.5 + 1))
print(f"\nTotal number of states: {final_answer}")
<<<6>>>