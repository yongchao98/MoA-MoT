import numpy as np

# Step 1 & 2: Define the Cartan matrices for the principal 2-blocks of A5 and C2.
# C(B0(A5)) from modular representation theory literature.
c_a5 = np.array([[2, 2],
                 [2, 4]])

# C(B0(C2)) is (2) as explained in the plan.
c_c2 = np.array([[2]])

# Step 3: The sum of entries for each matrix.
sum_a5 = np.sum(c_a5)
sum_c2 = np.sum(c_c2)

# The sum of entries of the tensor product is the product of the sums.
total_sum = sum_a5 * sum_c2

# Step 4: Display the logic and the final answer, showing the numbers used.
print(f"The Cartan matrix for the principal 2-block of A5 is:\n{c_a5}")
print(f"The sum of its entries is: {c_a5[0,0]} + {c_a5[0,1]} + {c_a5[1,0]} + {c_a5[1,1]} = {sum_a5}")
print("\n")
print(f"The Cartan matrix for the principal 2-block of C2 is:\n{c_c2}")
print(f"The sum of its entries is: {c_c2[0,0]}")
print("\n")
print("The sum of all entries in the Cartan matrix for the principal block of A5 x C2 is the product of these two sums.")
print(f"Final Calculation: ({c_a5[0,0]} + {c_a5[0,1]} + {c_a5[1,0]} + {c_a5[1,1]}) * {c_c2[0,0]} = {sum_a5} * {sum_c2} = {total_sum}")