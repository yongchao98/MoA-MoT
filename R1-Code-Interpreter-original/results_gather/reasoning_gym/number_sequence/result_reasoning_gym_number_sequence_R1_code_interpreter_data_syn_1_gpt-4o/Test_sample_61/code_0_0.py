# Define the sequence
sequence = [-5, -10, -20, -40, -80, -160, -320]

# Calculate the next term
last_term = sequence[-1]
last_difference = last_term - sequence[-2]
next_difference = last_difference * 2
next_term = last_term + next_difference

# Output the next term
print(next_term)