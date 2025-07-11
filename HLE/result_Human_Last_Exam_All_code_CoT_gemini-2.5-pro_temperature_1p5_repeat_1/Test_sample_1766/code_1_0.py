import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # For a p-electron, the orbital angular momentum quantum number is l=1.
    l = 1
    # For an electron, the spin quantum number is s=1/2.
    s = 0.5

    print(f"Solving for a p-electron with l={l} and s={s}.")
    print("-" * 60)

    # The total angular momentum quantum number j ranges from |l-s| to l+s.
    j_min = abs(l - s)
    j_max = l + s
    j_values = np.arange(j_min, j_max + 1, 1.0)

    # We will list the eigenvalues for each possible state |j, m_j>.
    # The eigenvalue of J^2 is j(j+1)hbar^2
    # The eigenvalue of J_z is m_j*hbar
    
    print("The common eigenvalues of J^2 and J_z are determined for each state |j, m_j>:\n")

    for j in j_values:
        # Calculate the eigenvalue for J^2 for the current j
        j2_eigenvalue = j * (j + 1)
        
        # The magnetic quantum number m_j ranges from -j to +j in steps of 1.
        mj_values = np.arange(-j, j + 1, 1.0)
        
        for mj in mj_values:
            print(f"For the state with quantum numbers j = {j} and m_j = {mj}:")
            # The eigenvalue of J^2 is printed showing the calculation
            print(f"  Eigenvalue of J^2 = j(j+1)hbar^2 = {j}({j}+1)hbar^2 = {j2_eigenvalue:.2f} hbar^2")
            # The eigenvalue of J_z is printed showing the calculation
            print(f"  Eigenvalue of J_z = m_j*hbar = {mj} hbar\n")

if __name__ == '__main__':
    solve_spin_orbit_coupling()
