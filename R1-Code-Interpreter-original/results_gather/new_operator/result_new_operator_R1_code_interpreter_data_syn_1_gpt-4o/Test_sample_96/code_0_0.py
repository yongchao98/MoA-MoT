# Step 1: Calculate 3@4
a = 3
b = 4
result_1 = 5 * a * b

# Step 2: Calculate (3@4)~5
c = 5
result_2 = 2 * result_1 + 3 * c

# Step 3: Calculate ((3@4)~5)#6
d = 6
final_result = 1 * result_2 * d

print(result_1, result_2, final_result)