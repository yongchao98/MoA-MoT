from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues for J^2 and J_z
    for a p-electron.
    """
    # 1. Define the quantum numbers for a p-electron.
    l = 1
    s = Fraction(1, 2)
    
    print(f"For a p-electron, the orbital quantum number l = {l}.")
    print(f"The spin quantum number s = {s}.")
    print("-" * 60)

    # 2. Determine the possible values for the total angular momentum quantum number, j.
    j_min = abs(l - s)
    j_max = l + s
    
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1
        
    print(f"The possible values for the total angular momentum quantum number j are {', '.join(map(str, j_values))}.")
    print("We will now find the eigenvalues for each possible j value.")
    print("Eigenvalue of J^2 = j(j+1)ħ²")
    print("Eigenvalue of J_z = m_j*ħ, where m_j = -j, -j+1, ..., +j")
    print("=" * 60)

    total_degeneracy = 0

    # 3. Loop through each j value to find the corresponding eigenvalues.
    for j in j_values:
        print(f"For the state with j = {j}:")
        
        # Calculate and display the J^2 eigenvalue, showing the formula.
        j_plus_1 = j + 1
        j2_eigenvalue = j * j_plus_1
        print(f"  The eigenvalue of J^2 is: {j} * ({j} + 1) ħ² = {j} * {j_plus_1} ħ² = {j2_eigenvalue} ħ²")

        # Determine the possible m_j values for this j.
        m_j_values = []
        current_m_j = -j
        while current_m_j <= j:
            m_j_values.append(current_m_j)
            current_m_j += 1
            
        total_degeneracy += len(m_j_values)
        
        # Display the corresponding J_z eigenvalues.
        print(f"  The corresponding eigenvalues of J_z are m_j * ħ, where m_j can be:")
        m_j_str = ", ".join(map(str, m_j_values))
        print(f"  m_j = {m_j_str}")
        print("-" * 60)
        
    print(f"The total degeneracy (total number of states) is the sum of the multiplicities (2j+1) for each j.")
    print(f"Initial state degeneracy: (2l+1)*(2s+1) = (2*{l}+1)*(2*{s}+1) = {int((2*l+1)*(2*s+1))}")
    print(f"Final state degeneracy: sum over j of (2j+1) = {total_degeneracy}")
    print("The degeneracy is conserved, as expected.")


solve_spin_orbit_coupling()

<<<6>>>