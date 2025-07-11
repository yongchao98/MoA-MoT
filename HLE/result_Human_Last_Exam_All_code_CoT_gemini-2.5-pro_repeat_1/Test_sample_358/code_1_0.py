import numpy as np

# Step 1: Define the decomposition matrix for the principal 2-block of A5.
# The block has 4 ordinary characters and 2 simple modules.
D_A5 = np.array([
    [1, 0],
    [1, 1],
    [1, 1],
    [1, 2]
])

# Step 2: Calculate the Cartan matrix for A5 using C = D^T * D.
C_A5 = D_A5.T @ D_A5

# Step 3: Calculate the sum of all entries in the Cartan matrix for A5.
sum_A5 = np.sum(C_A5)

print("For the principal 2-block of A5:")
print("Decomposition matrix D_A5:")
print(D_A5)
print("\nCartan matrix C_A5 = D_A5^T * D_A5:")
print(C_A5)
print(f"\nSum of entries in C_A5 = {C_A5[0,0]} + {C_A5[0,1]} + {C_A5[1,0]} + {C_A5[1,1]} = {int(sum_A5)}")

print("-" * 30)

# Step 4: Define the Cartan matrix for the principal 2-block of C2.
# It's a 1x1 matrix (2).
C_C2 = np.array([[2]])
sum_C2 = np.sum(C_C2)

print("For the principal 2-block of C2:")
print("Cartan matrix C_C2:")
print(C_C2)
print(f"\nSum of entries in C_C2 = {int(sum_C2)}")

print("-" * 30)

# Step 5: The total sum is the product of the individual sums.
total_sum = sum_A5 * sum_C2

print("The sum of all entries in the Cartan matrix of the principal block of A5 x C2")
print(f"is the product of the individual sums:")
print(f"{int(sum_A5)} * {int(sum_C2)} = {int(total_sum)}")

# Final answer in the required format
# print(f"\n<<<{int(total_sum)}>>>")