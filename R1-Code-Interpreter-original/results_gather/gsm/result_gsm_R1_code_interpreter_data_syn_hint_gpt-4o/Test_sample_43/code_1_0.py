# Initial popping rate
initial_rate = 360935

# Calculate the number of kernels popped in each interval
first_interval = initial_rate
second_interval = 3 * initial_rate
third_interval = 4 * initial_rate
fourth_interval = 2 * initial_rate
residual_popping = round(0.5 * initial_rate)  # Round to the nearest whole number

# Total number of popped kernels
total_popped = first_interval + second_interval + third_interval + fourth_interval + residual_popping

# Output the total number of popped kernels
print(total_popped)