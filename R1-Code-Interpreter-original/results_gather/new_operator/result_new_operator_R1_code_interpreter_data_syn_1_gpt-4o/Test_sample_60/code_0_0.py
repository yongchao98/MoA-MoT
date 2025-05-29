# Step 1: Calculate 3~4
a = 3
b = 4
result_1 = 5 * a + 3 * b

# Step 2: Calculate (3~4)%5
a = result_1
b = 5
result_2 = 3 * a + 2 * b

# Step 3: Calculate (3~4)%5$6
a = result_2
b = 6
result_3 = (a + b) / b

print(result_1, result_2, result_3)