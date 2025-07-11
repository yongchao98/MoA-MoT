import numpy as np

# Step 1: Define the Cartan matrix for the principal 2-block of A5.
# This matrix is a known result from modular representation theory.
cartan_matrix_A5 = np.array([
    [4, 2, 2],
    [2, 2, 1],
    [2, 1, 2]
])

# Step 2: Calculate the sum of all entries in the Cartan matrix for A5.
sum_A5 = np.sum(cartan_matrix_A5)

# Step 3: Define the sum of entries for the Cartan matrix of the principal 2-block of C2.
# For C2 in characteristic 2, the Cartan matrix is (2), so the sum is 2.
sum_C2 = 2

# Step 4: The total sum is the product of the individual sums.
total_sum = sum_A5 * sum_C2

# Step 5: Print the final equation as requested.
# The output shows the sum for A5, the sum for C2, and the final product.
print(f"The sum of entries for the A5 principal block's Cartan matrix is {int(sum_A5)}.")
print(f"The sum of entries for the C2 principal block's Cartan matrix is {sum_C2}.")
print("The total sum for the principal block of A5 x C2 is the product of these two sums:")
print(f"{int(sum_A5)} * {sum_C2} = {int(total_sum)}")
