# Step 1: Calculate 3~4
a = 3
b = 4
result_1 = 5 * a + 4 * b

# Step 2: Calculate (3~4)%5
b = 5
result_2 = 2 * result_1 * b

# Step 3: Calculate ((3~4)%5)@6
b = 6
result_3 = 4 * result_2 + 4 * b

print(result_1, result_2, result_3)