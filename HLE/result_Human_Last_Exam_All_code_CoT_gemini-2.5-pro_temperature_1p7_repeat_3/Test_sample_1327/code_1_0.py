# The sequence of numbers
n = [2, 11, 23, 51, 119]

# Step 1: Establish the primary pattern: n_i+1 = 2 * n_i + a_i
# We calculate the sequence 'a'
a = []
for i in range(len(n) - 1):
    a_i = n[i+1] - 2 * n[i]
    a.append(a_i)
    # print(f"{n[i+1]} = 2 * {n[i]} + {a_i}")

# The calculated sequence 'a' is [7, 1, 5, 17]

# Step 2: Find the pattern in the sub-sequence of 'a' starting from the second term: [1, 5, 17]
# We calculate the differences between its terms
diff1 = a[2] - a[1]  # 5 - 1
diff2 = a[3] - a[2]  # 17 - 5
# print(f"The differences are {diff1} and {diff2}.")

# Step 3: Identify the pattern in the differences [4, 12]
# It's a geometric progression with a ratio of 3.
ratio = diff2 // diff1
# print(f"The ratio is {ratio}.")

# Step 4: Predict the next difference and the next term in 'a'
next_diff = diff2 * ratio
next_a = a[-1] + next_diff
# print(f"The next term in sequence 'a' is {a[-1]} + {next_diff} = {next_a}")

# Step 5: Calculate the final missing number in the original sequence
last_known_n = n[-1]
final_result = 2 * last_known_n + next_a

# Print the final calculation and the answer
print(f"The pattern is that each number is twice the previous number plus an increasing value.")
print(f"The values being added are determined by a sub-pattern.")
print(f"The next value to be added is {next_a}.")
print(f"Therefore, the final equation is:")
print(f"2 * {last_known_n} + {next_a} = {final_result}")
<<<A>>>