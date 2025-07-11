import numpy as np

def solve_spin_orbit_eigenvalues():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    We will use atomic units where ħ (h-bar) is 1, but we will show it in the formulas for clarity.
    """
    # Step 1: Define the quantum numbers for a p-electron.
    l = 1.0  # Orbital angular momentum quantum number for a p-orbital
    s = 0.5  # Spin quantum number for an electron

    print(f"For a p-electron, the orbital angular momentum quantum number is l = {int(l)}.")
    print(f"The spin quantum number is s = {s}.\n")

    # Step 2: Determine the possible values for the total angular momentum quantum number, j.
    # j ranges from |l-s| to l+s in integer steps.
    j_min = abs(l - s)
    j_max = l + s
    j_values = np.arange(j_min, j_max + 1, 1.0)

    print(f"The possible total angular momentum quantum numbers (j) are: {', '.join(map(str, j_values))}")
    print("--------------------------------------------------\n")
    
    final_answer_parts = []

    # Step 3 & 4: For each j, find m_j and calculate the eigenvalues.
    for j in j_values:
        print(f"Analyzing the state with j = {j}:")
        
        # Calculate and display the eigenvalue for J^2 = j(j+1)ħ²
        j_squared_eigenvalue = j * (j + 1)
        print(f"  - The eigenvalue of the J^2 operator is j(j+1)ħ².")
        print(f"    Calculation: {j} * ({j} + 1)ħ² = {j_squared_eigenvalue}ħ²\n")

        # Determine the possible m_j values, which range from -j to j in integer steps.
        m_j_values = np.arange(-j, j + 1, 1.0)
        
        # Display the eigenvalues for J_z = m_jħ
        print(f"  - For this j, the possible m_j quantum numbers are: {', '.join(map(str, m_j_values))}.")
        print(f"  - The corresponding eigenvalues of the J_z operator are m_jħ.")
        
        print("--------------------------------------------------\n")

        # Build the final summary answer string
        part = (f"For j={j}, the J^2 eigenvalue is {j_squared_eigenvalue}ħ^2 and the J_z eigenvalues are "
                f"{{{', '.join(map(str, m_j_values))}}}ħ.")
        final_answer_parts.append(part)

    # Combine parts for the final answer
    final_answer = " ".join(final_answer_parts)
    print(f"<<<{final_answer}>>>")


# Run the function
solve_spin_orbit_eigenvalues()