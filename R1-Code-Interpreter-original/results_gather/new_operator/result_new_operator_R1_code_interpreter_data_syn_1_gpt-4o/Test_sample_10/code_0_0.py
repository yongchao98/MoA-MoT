# Step 1: Calculate 3~4
a = 3
b = 4
result_3_tilde_4 = (a + b) / b

# Step 2: Calculate (3~4)#5
c = 5
result_3_tilde_4_hash_5 = 5 * result_3_tilde_4 + 5 * c

# Step 3: Calculate (3~4)#5@6
d = 6
result_final = result_3_tilde_4_hash_5 * d

print(result_3_tilde_4, result_3_tilde_4_hash_5, result_final)