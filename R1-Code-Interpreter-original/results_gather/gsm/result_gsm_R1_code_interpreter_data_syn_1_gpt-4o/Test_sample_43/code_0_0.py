# Initial popping rate
initial_rate = 360935

# Calculate the number of kernels popped in each interval
first_30_seconds = initial_rate
second_30_seconds = 3 * initial_rate
third_30_seconds = 4 * initial_rate
fourth_30_seconds = (1/2) * third_30_seconds
residual_popping = (1/4) * fourth_30_seconds

# Total number of popped kernels
total_popped = first_30_seconds + second_30_seconds + third_30_seconds + fourth_30_seconds + residual_popping

# Print the total number of popped kernels
print(total_popped)