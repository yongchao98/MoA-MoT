# The given sequence
sequence = [2, 11, 23, 51, 119]

# The first term used in the pattern calculation
a2 = sequence[1] 
# The second term
a3 = sequence[2] 
# The third term
a4 = sequence[3] 
# The last known term
a5 = sequence[4] 

# Calculate the sequence of subtracted numbers (c_n)
c3 = 3 * a2 - a3
c4 = 3 * a3 - a4
c5 = 3 * a4 - a5

# Calculate the differences in the c_n sequence
diff1 = c4 - c3
diff2 = c5 - c4

# Predict the next difference (doubles each time)
next_diff = diff2 * 2

# Predict the next subtracted number
next_c = c5 + next_diff

# Calculate the next term in the main sequence
next_term = 3 * a5 - next_c

# Output the explanation and the final equation
print("The pattern is a(n) = 3 * a(n-1) - c(n), where c(n) follows its own pattern.")
print(f"The sequence of subtracted numbers starts with: {c3}, {c4}, {c5}, ...")
print(f"The differences in this sequence are {diff1}, {diff2}, ... which double each time.")
print(f"The next difference is {diff2} * 2 = {next_diff}.")
print(f"So the next number to subtract is {c5} + {next_diff} = {next_c}.")
print("\nThe final calculation is:")
print(f"{3} * {a5} - {next_c} = {next_term}")