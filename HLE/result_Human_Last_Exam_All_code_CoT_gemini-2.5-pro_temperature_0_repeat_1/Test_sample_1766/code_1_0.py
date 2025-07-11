from fractions import Fraction
import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    
    A p-electron has orbital angular momentum l=1 and spin s=1/2.
    The total angular momentum J is the sum of L and S.
    The eigenvalues of J^2 are j(j+1) h_bar^2.
    The eigenvalues of J_z are m_j h_bar.
    """
    # Define quantum numbers for a p-electron
    l = 1
    s = Fraction(1, 2)

    print(f"Solving for a p-electron with l={l} and s={s}.\n")

    # Calculate possible j values from |l-s| to l+s
    j_min = abs(l - s)
    j_max = l + s
    j_values = [Fraction(j) for j in np.arange(float(j_min), float(j_max) + 1, 1.0)]

    # Iterate through each possible j value
    for j in sorted(j_values, reverse=True):
        print(f"--- For states with total angular momentum j = {j} ---")
        
        # Calculate the eigenvalue for J^2
        j2_eigenvalue = j * (j + 1)
        
        # m_j values range from -j to j in integer steps
        m_j_values = [Fraction(m) for m in np.arange(float(-j), float(j) + 1, 1.0)]
        
        # Iterate through each m_j for the given j
        for m_j in m_j_values:
            state_ket = f"|j={j}, m_j={m_j}>"
            
            # Print the eigenvalue equation for J^2
            print(f"For state {state_ket}:")
            print(f"  J^2 {state_ket} = {j}*({j}+1) h_bar^2 {state_ket} = {j2_eigenvalue} h_bar^2 {state_ket}")
            
            # Print the eigenvalue equation for J_z
            print(f"  J_z {state_ket} = {m_j} h_bar {state_ket}\n")

solve_spin_orbital_coupling()