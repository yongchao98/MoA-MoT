import fractions

def solve_spin_orbital_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z for a p-electron.
    The eigenvalues of J^2 are hbar^2 * j*(j+1) and for J_z are hbar * m_j.
    This script calculates the numerical factors j*(j+1) and m_j.
    """
    # For a p-electron, orbital angular momentum quantum number l = 1
    l = 1
    # For an electron, spin quantum number s = 1/2
    s = 0.5

    print(f"Solving for a p-electron with l={l} and s={s}:\n")

    # The total angular momentum quantum number j can take values from |l-s| to l+s
    j_min = abs(l - s)
    j_max = l + s

    # Generate the list of possible j values
    j_values = []
    current_j = j_min
    while current_j <= j_max:
        j_values.append(current_j)
        current_j += 1

    # Store results for the final answer
    eigenvalue_pairs = []

    # Iterate through each possible j manifold
    for j in j_values:
        j_frac = fractions.Fraction(j).limit_denominator()
        print(f"--- For the j = {j_frac} manifold ---")

        # The eigenvalue of J^2 is hbar^2 * j * (j + 1)
        j_squared_eigenvalue_factor = j * (j + 1)
        j2_frac = fractions.Fraction(j_squared_eigenvalue_factor).limit_denominator()

        # The magnetic quantum number m_j ranges from -j to +j in integer steps
        m_j_values = []
        current_m_j = -j
        while current_m_j <= j:
            m_j_values.append(current_m_j)
            current_m_j += 1

        print(f"The J^2 eigenvalue is constant for this manifold: {j2_frac} * hbar^2")
        print("The corresponding J_z eigenvalues are:")
        
        # For each m_j, the J_z eigenvalue is hbar * m_j
        for m_j in m_j_values:
            mj_frac = fractions.Fraction(m_j).limit_denominator()
            
            # Add to result list
            eigenvalue_pairs.append(((str(j2_frac)), str(mj_frac)))

            print(f"  For the state |j={j_frac}, m_j={mj_frac}>:")
            # Output each number in the final equation per user request
            print(f"    Eigenvalue of J^2 = ({j_frac}) * ({j_frac} + 1) * hbar^2 = {j2_frac} * hbar^2")
            print(f"    Eigenvalue of J_z = {mj_frac} * hbar")
        print() # Add a newline for better separation
    
    # This part is for the final answer block and will not be printed to console
    final_answer_str = "The pairs of eigenvalues (for J^2, for J_z) in units of hbar^2 and hbar respectively are: "
    final_answer_str += ", ".join([f"({p[0]}, {p[1]})" for p in eigenvalue_pairs])
    return final_answer_str

# Execute the function and print the results
final_answer = solve_spin_orbital_coupling()
# The final answer is wrapped for the platform, but the script itself produces the detailed output.
# print(f"\n<<<final_answer>>>")
