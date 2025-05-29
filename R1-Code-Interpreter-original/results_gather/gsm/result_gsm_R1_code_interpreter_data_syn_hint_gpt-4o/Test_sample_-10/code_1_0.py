# Initial popping rate
initial_rate = 360935

# Kernels popped in each interval
first_interval = initial_rate
second_interval = 3 * initial_rate
third_interval = 4 * initial_rate
fourth_interval = 2 * initial_rate
residual_popping = round(0.5 * initial_rate)

# Total kernels popped
total_kernels = first_interval + second_interval + third_interval + fourth_interval + residual_popping

print(total_kernels)