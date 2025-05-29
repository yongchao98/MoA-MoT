# Step 1: Calculate 3@4
a = 3
b = 4
result_3_at_4 = (a + b) / b

# Step 2: Calculate (3@4)#5
a = result_3_at_4
b = 5
result_3_at_4_hash_5 = 2 * a + 3 * b

# Step 3: Calculate ((3@4)#5)$6
a = result_3_at_4_hash_5
b = 6
final_result = (a + b) / b

print(final_result)