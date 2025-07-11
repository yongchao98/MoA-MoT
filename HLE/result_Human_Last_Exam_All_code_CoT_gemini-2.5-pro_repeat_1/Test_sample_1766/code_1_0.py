from fractions import Fraction

# For a p-electron, the orbital angular momentum quantum number l = 1.
# The spin angular momentum quantum number for an electron is s = 1/2.
l = 1
s = Fraction(1, 2)

print(f"Solving for common eigenvalues of J^2 and J_z for a p-electron (l={l}, s={s}).")
print("Eigenvalues are given in units of ħ=1.\n")

# The total angular momentum quantum number j can take values from |l-s| to l+s.
j_min = abs(l - s)
j_max = l + s

j_values = []
current_j = j_min
while current_j <= j_max:
    j_values.append(current_j)
    current_j += 1

print(f"The possible total angular momentum quantum numbers (j) are: {', '.join(map(str, j_values))}\n")
print("The common eigenvalues are found for each common eigenstate |j, m_j>.")
print("-" * 60)

# Store the eigenvalue pairs for the final answer
eigenvalue_pairs = []

# For each j, find the eigenvalues of J^2 and J_z.
for j in j_values:
    # The eigenvalue of J^2 is j*(j+1)
    j_squared_eigenvalue = j * (j + 1)
    
    # The magnetic quantum number m_j can take values from -j to +j in steps of 1.
    m_j_values = []
    current_m_j = -j
    while current_m_j <= j:
        m_j_values.append(current_m_j)
        current_m_j += 1
    
    print(f"For the subspace with j = {j}:")
    
    # The eigenvalue of J^2 is the same for all states with the same j.
    # The prompt requires showing the numbers in the equation.
    print(f"  The eigenvalue of J^2 is j(j+1) = {j} * ({j} + 1) = {j_squared_eigenvalue}")
    
    print(f"  The possible m_j values are: {', '.join(map(str, m_j_values))}")
    print("  The corresponding common eigenvalue pairs (J^2, J_z) are:")
    
    for m_j in m_j_values:
        # The eigenvalue of J_z is m_j
        j_z_eigenvalue = m_j
        print(f"    For state |j={j}, m_j={m_j}>:\t({j_squared_eigenvalue}, {j_z_eigenvalue})")
        eigenvalue_pairs.append((j_squared_eigenvalue, j_z_eigenvalue))
        
    print("-" * 60)