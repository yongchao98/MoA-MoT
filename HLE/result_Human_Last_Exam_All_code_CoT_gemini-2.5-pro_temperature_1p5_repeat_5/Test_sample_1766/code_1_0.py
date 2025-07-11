import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    The eigenvalues are the physical observables for a system in a state defined
    by the quantum numbers j and m_j.
    """
    # For a p-electron, the orbital angular momentum quantum number is 1.
    l = 1
    # For an electron, the spin quantum number is 1/2.
    s = 0.5

    print(f"Solving for a p-electron with l = {l} and s = {s}.")
    print("The total angular momentum is J = L + S.")
    print("We need to find the common eigenvalues of the operators J^2 and J_z.")
    print("-" * 60)

    # The rules of angular momentum addition (Clebsch-Gordan series) state that
    # the total angular momentum quantum number j can take values from |l-s| to l+s.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print(f"The possible values for the total angular momentum quantum number j are: {list(j_values)}\n")

    # Iterate through each possible value of j to find the corresponding eigenvalues.
    for j in j_values:
        # The eigenvalue of the J^2 operator is j(j+1)hbar^2.
        j_squared_eigenvalue = j * (j + 1)
        
        print(f"--- For the state with j = {j} ---")
        print(f"The eigenvalue of J^2 is given by j(j+1)ħ².")
        # We explicitly show the calculation for clarity.
        print(f"  Eigenvalue of J^2 = {j} * ({j} + 1) * ħ² = {j_squared_eigenvalue:.2f} * ħ²")

        # For a given j, the magnetic quantum number m_j can take values from -j to +j.
        mj_values = np.arange(-j, j + 1, 1)
        
        print(f"The possible values for the magnetic quantum number m_j are: {list(mj_values)}")
        print("The corresponding J_z eigenvalues (m_j * ħ) for each common eigenstate |j, m_j> are:")

        # For each state |j, m_j>, there is a specific eigenvalue for J_z.
        for mj in mj_values:
            # The eigenvalue of the J_z operator is m_j * hbar.
            print(f"  For state |j={j}, m_j={mj:+.1f}>, the eigenvalue of J_z is {mj:+.1f} * ħ")
        
        print("-" * 60)

if __name__ == "__main__":
    solve_spin_orbit_coupling()