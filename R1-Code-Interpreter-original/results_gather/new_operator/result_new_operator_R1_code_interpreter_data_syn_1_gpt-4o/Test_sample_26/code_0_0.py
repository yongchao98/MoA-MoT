# Step 1: Calculate 3#4
a = 3
b = 4
result_3_hash_4 = (a + b) / b

# Step 2: Calculate (3#4)@5
a = result_3_hash_4
b = 5
result_3_hash_4_at_5 = 5 * a + 5 * b

# Step 3: Calculate (3#4)@5~6
a = result_3_hash_4_at_5
b = 6
result_final = 3 * a * b

print(result_3_hash_4, result_3_hash_4_at_5, result_final)