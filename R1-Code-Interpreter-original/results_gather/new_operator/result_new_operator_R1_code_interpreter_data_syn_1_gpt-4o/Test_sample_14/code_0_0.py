# Step 1: Calculate 3$4
a = 3
b = 4
result_1 = 2 * a + 5 * b

# Step 2: Calculate (3$4)#5
b = 5
result_2 = 5 * result_1 * b

# Step 3: Calculate ((3$4)#5)~6
b = 6
result_3 = (result_2 + b) / b

print(result_1, result_2, result_3)