import numpy as np

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    """
    # For a p-electron, orbital angular momentum quantum number l=1
    l = 1
    # For an electron, spin angular momentum quantum number s=1/2
    s = 0.5

    print(f"For a p-electron with l = {l} and s = {s}:\n")

    # The total angular momentum quantum number j can take values from |l-s| to l+s.
    j_values = np.arange(abs(l - s), l + s + 1, 1)

    print("The eigenvalues of J^2 are j(j+1)hbar^2 and for J_z are m_j*hbar.")
    print("Setting hbar = 1, the common eigenvalues are:\n")

    # Store the results for the final answer
    eigenvalue_pairs = []

    for j in j_values:
        # Eigenvalue of J^2 is j*(j+1)
        j2_eigenvalue = j * (j + 1)
        
        print("-----------------------------------------")
        print(f"For the total angular momentum state j = {j}:")
        
        # Print the J^2 eigenvalue calculation
        print(f"The eigenvalue of J^2 is {j} * ({j} + 1) = {j2_eigenvalue}")
        
        # m_j ranges from -j to +j in steps of 1
        m_j_values = np.arange(-j, j + 1, 1)
        
        print(f"The corresponding eigenvalues of J_z (m_j) are: {list(m_j_values)}")
        print("-----------------------------------------\n")
        
        for m_j in m_j_values:
            eigenvalue_pairs.append((j2_eigenvalue, m_j))

    # This part is just to generate the final answer string, not for the user to see in the script's output
    final_answer = ", ".join([f"({p[0]}, {p[1]})" for p in eigenvalue_pairs])
    # print(f"<<<{final_answer}>>>")


solve_spin_orbital_coupling()