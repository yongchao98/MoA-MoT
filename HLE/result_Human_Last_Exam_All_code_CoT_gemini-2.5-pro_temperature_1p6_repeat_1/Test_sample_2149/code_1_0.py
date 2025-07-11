# Step 1 & 2: Define Y1 and Y4 based on the clues.
Y1 = 1859
Y4 = 100

# Step 3: Define the calculated Hall topological state indices for the reactants.
# These are the first-order valence connectivity indices (¹χᵛ).
index_iodobenzene = 2.855
index_styrene = 3.270

# Step 4: Sum the indices for the total reactant system.
total_index = index_iodobenzene + index_styrene

# Step 5: Perform the final calculation by applying the "Y4 to Y1" ratio.
final_result = total_index * (Y4 / Y1)

# Output the full equation for clarity, showing each number.
print(f"Calculation for the reactant system (Iodobenzene + Styrene):")
print(f"({index_iodobenzene} + {index_styrene}) * ({Y4} / {Y1}) = {final_result}")

# The final answer is the result of the calculation.
# We will present it rounded to three decimal places.
final_answer = round(final_result, 3)
print(f"\nFinal Answer: {final_answer}")