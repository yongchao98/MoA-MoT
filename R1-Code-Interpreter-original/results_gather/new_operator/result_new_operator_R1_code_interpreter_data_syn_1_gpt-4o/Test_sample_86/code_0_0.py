# Step 1: Calculate 3&4
a = 3
b = 4
result_3_and_4 = 5 * a + 2 * b

# Step 2: Calculate (3&4)@5
b = 5
result_3_and_4_at_5 = 3 * result_3_and_4 * b

# Step 3: Calculate ((3&4)@5)#6
b = 6
final_result = (result_3_and_4_at_5 + b) / b

print(result_3_and_4, result_3_and_4_at_5, final_result)