import numpy as np
from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # For a p-electron, the orbital angular momentum quantum number is l=1.
    l = 1
    # For an electron, the spin quantum number is s=1/2.
    s = 0.5

    print(f"Solving for a p-electron with l={l} and s={Fraction(s)}.\n")

    # The rules of angular momentum addition state that the total angular
    # momentum quantum number j can take values from |l-s| to l+s.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print("The eigenvalue equations are of the form:")
    print("J^2 |j, m_j> = j(j+1) * h_bar^2 |j, m_j>")
    print("J_z |j, m_j> = m_j * h_bar |j, m_j>\n")

    # Iterate through each possible value of j
    for j in j_values:
        j_frac = Fraction(j).limit_denominator()
        
        # Calculate the eigenvalue for J^2, which depends only on j.
        j2_eigenvalue = Fraction(j * (j + 1)).limit_denominator()
        
        print(f"--- For states with total angular momentum j = {j_frac} ---")
        
        # For a given j, the magnetic quantum number m_j runs from -j to +j.
        mj_values = np.arange(-j, j + 1, 1)
        
        # Iterate through each possible value of m_j for the given j
        for mj in mj_values:
            mj_frac = Fraction(mj).limit_denominator()
            
            # Print the state and its corresponding eigenvalues
            print(f"For the state |j={j_frac}, m_j={mj_frac}>:")
            print(f"  J^2 |{j_frac}, {mj_frac}> = {j2_eigenvalue} * h_bar^2 |{j_frac}, {mj_frac}>")
            print(f"  J_z |{j_frac}, {mj_frac}> = {mj_frac} * h_bar |{j_frac}, {mj_frac}>\n")

if __name__ == "__main__":
    solve_spin_orbit_coupling()