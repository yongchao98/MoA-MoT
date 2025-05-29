# Step 1: Calculate 3#4
a = 3
b = 4
result_1 = 4 * a * b

# Step 2: Calculate (3#4)~5
c = 5
result_2 = 3 * result_1 + 5 * c

# Step 3: Calculate ((3#4)~5)%6
d = 6
result_3 = 4 * result_2 * d

print(result_1, result_2, result_3)