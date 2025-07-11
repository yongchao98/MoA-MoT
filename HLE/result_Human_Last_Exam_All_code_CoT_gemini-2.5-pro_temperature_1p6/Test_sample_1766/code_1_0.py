from fractions import Fraction

# Step 1: Explain the setup
print("Solving for the common eigenvalues of J^2 and J_z for a p-electron.")
print("A p-electron has orbital angular momentum quantum number l=1.")
print("An electron has spin angular momentum quantum number s=1/2.")
print("We will work in units where ħ = 1.")
print("-" * 60)

# Step 2: Define quantum numbers
l = 1
s = 0.5

# Step 3: Determine the possible values for the total angular momentum quantum number, j
# j ranges from |l-s| to l+s in integer steps.
j_values = [abs(l - s), l + s]

print(f"The possible values for the total angular momentum quantum number j are: {', '.join(map(str, j_values))}\n")

# Step 4: Loop through each possible j value to find the eigenvalues
for j in j_values:
    # Eigenvalue of J^2 is j(j+1)
    j2_eigenvalue = j * (j + 1)
    
    # Possible values for m_j, the z-component of J, range from -j to j in steps of 1
    # Generate m_j values from -j to j
    num_mj_steps = int(2 * j + 1)
    mj_values = [-j + i for i in range(num_mj_steps)]
    
    j_frac = Fraction(j)
    j2_eigenvalue_frac = Fraction(j2_eigenvalue)

    print(f"For the case j = {j_frac}:")
    print(f"The eigenvalue of J^2 is calculated as j*(j+1):")
    print(f"  J^2 eigenvalue = {j_frac} * ({j_frac} + 1) = {j2_eigenvalue_frac}")
    
    print("\nThe corresponding common eigenstates |j, m_j> have the following J^2 and J_z eigenvalues:")
    
    # Loop through all possible m_j values for the given j
    for mj in mj_values:
        # Eigenvalue of J_z is m_j
        jz_eigenvalue = mj
        
        mj_frac = Fraction(mj)
        jz_eigenvalue_frac = Fraction(jz_eigenvalue)
        
        print(f"  For state |j={j_frac}, m_j={mj_frac}>:")
        # Final eigenvalue equations
        print(f"    J^2 |{j_frac}, {mj_frac}> = {j2_eigenvalue_frac} |{j_frac}, {mj_frac}>")
        print(f"    J_z |{j_frac}, {mj_frac}> = {jz_eigenvalue_frac} |{j_frac}, {mj_frac}>")
        
    print("-" * 60)