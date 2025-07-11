# Plan:
# 1. Identify the sexes for the insects in pair A based on abdominal shape.
#    - Left (pointed abdomen) is female.
#    - Right (blunter, multi-spined tip) is male.
#    - Pair A is F, M, which corresponds to index 4.
# 2. Identify the sexes for the insects in pair B based on antennae shape.
#    - Left (curled antennae) is male.
#    - Right (straight antennae) is female.
#    - Pair B is M, F, which corresponds to index 3.
# 3. Identify the sexes for the insects in pair C based on antennae length.
#    - Left (very long antennae) is male.
#    - Right (short antennae) is female.
#    - Pair C is M, F, which corresponds to index 3.
# 4. Combine the indices into the required format "A,B,C".

# Final indices for pairs A, B, and C
index_A = 4
index_B = 3
index_C = 3

# Print the final answer in the specified format.
print(f"{index_A},{index_B},{index_C}")
