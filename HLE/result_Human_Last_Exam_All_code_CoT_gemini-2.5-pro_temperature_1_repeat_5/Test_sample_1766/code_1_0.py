from fractions import Fraction

# Step 1: Explain the problem setup
print("This script solves for the common eigenvalues of J^2 and J_z for a p-electron.")
print("A p-electron has orbital angular momentum quantum number l=1.")
print("An electron has spin quantum number s=1/2.")
print("The total angular momentum J is the sum of orbital (L) and spin (S) angular momentum: J = L + S.")
print("Since J^2 and J_z commute, we can find common eigenstates |j, m_j>.")
print("-" * 70)

# Step 2: Define quantum numbers using the Fraction class for pretty printing
l = 1
s = Fraction(1, 2)

# Step 3: Calculate the possible values for the total angular momentum quantum number j
# j ranges from |l-s| to l+s in integer steps.
j_min = abs(l - s)
j_max = l + s
# In this specific case, the steps are of size 1, so we get two values.
j_values = [j_min, j_max]

print(f"Possible values for the total angular momentum quantum number j are: {[str(j) for j in j_values]}.\n")

# Step 4: Iterate through all possible states |j, m_j> and print the eigenvalue equations
print("The eigenvalue equations for each common eigenstate |j, m_j> are:")
for j in j_values:
    # Calculate the eigenvalue for J^2, which is j*(j+1)
    j2_eigenvalue = j * (j + 1)
    
    # Calculate the possible values for m_j, which range from -j to +j in steps of 1
    current_mj = -j
    mj_values = []
    while current_mj <= j:
        mj_values.append(current_mj)
        current_mj += 1
        
    # Print the results for the current j-multiplet
    for mj in mj_values:
        # The eigenvalue for J_z is simply m_j
        jz_eigenvalue = mj
        
        # Print the final equations with all numbers substituted
        print(f"\nFor the state |j={j}, m_j={mj}>:")
        print(f"  J^2 |{j}, {mj}> = {j}({j}+1)ħ^2 |{j}, {mj}> = {j2_eigenvalue}ħ^2 |{j}, {mj}>")
        print(f"  J_z |{j}, {mj}> = {jz_eigenvalue}ħ |{j}, {mj}>")