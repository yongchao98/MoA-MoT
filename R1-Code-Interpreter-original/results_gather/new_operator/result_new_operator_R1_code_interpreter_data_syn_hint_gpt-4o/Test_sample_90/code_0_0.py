# Step 1: Calculate 3%4
a = 3
b = 4
result_1 = (a + b) / b

# Step 2: Calculate (3%4)$5
c = 5
result_2 = (result_1 + c) / c

# Step 3: Calculate (3%4)$5~6
d = 6
result_3 = 1 * result_2 + 5 * d

print(result_1, result_2, result_3)