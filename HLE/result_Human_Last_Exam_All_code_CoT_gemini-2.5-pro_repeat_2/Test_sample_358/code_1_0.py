import numpy as np

# Step 1: Explain the theoretical background.
print("The group is G = A5 x C2, and the characteristic of the field k is 2.")
print("The sum of entries of the Cartan matrix for the principal block of kG")
print("is the product of the sums of entries for the principal blocks of kA5 and kC2.")
print("-" * 60)

# Step 2: Analyze the C2 case.
print("For the principal block of kC2 in characteristic 2:")
# In characteristic 2, the group algebra kC2 has one simple module. Its projective
# cover has this simple module as a composition factor twice.
# Thus, the Cartan matrix is C_C2 = [[2]].
sum_c2 = 2
print("The Cartan matrix is a 1x1 matrix with the value 2.")
print(f"The sum of its entries is: {sum_c2}")
print("-" * 60)

# Step 3: Analyze the A5 case.
print("For the principal block of kA5 in characteristic 2:")
# This part uses known results from the modular representation theory of A5.
# The principal block of kA5 contains 3 simple modules. The ordinary characters
# in the block (with dimensions 1, 3, 3, 5) decompose into these simple modules.
# This is described by the decomposition matrix D_A5.
# Rows correspond to ordinary characters, columns to simple modules.
D_A5 = np.array([
    [1, 0, 0],
    [1, 1, 0],
    [1, 0, 1],
    [1, 1, 1]
])

print("The decomposition matrix D_A5 for the principal block is:")
print(D_A5)

# The Cartan matrix C is calculated as C = D^T * D.
C_A5 = D_A5.T @ D_A5
print("\nThe Cartan matrix C_A5 = D_A5^T * D_A5 is:")
print(C_A5)

# Sum the entries of the calculated Cartan matrix C_A5.
c_a5_entries = C_A5.flatten()
sum_a5 = int(np.sum(c_a5_entries))

# Display the summation equation for the entries of C_A5.
equation_str = " + ".join(map(str, c_a5_entries))
print("\nThe sum of the entries of the Cartan matrix is:")
print(f"{equation_str} = {sum_a5}")
print("-" * 60)

# Step 4: Calculate the final result.
total_sum = sum_a5 * sum_c2
print("The total sum of entries for the Cartan matrix of k[A5 x C2] is the product of the two sums:")
print(f"Sum(C_A5) * Sum(C_C2) = {sum_a5} * {sum_c2} = {total_sum}")

print(f"\n<<< {total_sum} >>>")