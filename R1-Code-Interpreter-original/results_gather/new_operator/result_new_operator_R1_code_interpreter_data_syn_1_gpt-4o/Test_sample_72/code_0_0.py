# Step 1: Calculate 3&4
a = 3
b = 4
c = 5
d = 6

result_3_and_4 = 3 * a + 2 * b

# Step 2: Calculate (3&4)~5
result_3_and_4_tilde_5 = 4 * result_3_and_4 * c

# Step 3: Calculate ((3&4)~5)@6
final_result = (result_3_and_4_tilde_5 + d) / d

print(result_3_and_4, result_3_and_4_tilde_5, final_result)