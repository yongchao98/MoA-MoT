# Step 1: Define the parameters for the calculation.
# The chosen test is March RAW, which has a complexity of 13N.
complexity_factor = 13
# RAM size is 1Mbit.
ram_size_bits = 1_000_000
# Cycle time is 5 nanoseconds.
cycle_time_ns = 5

# Step 2: Calculate the total number of operations.
total_operations = complexity_factor * ram_size_bits

# Step 3: Calculate the total test time in nanoseconds.
total_time_ns = total_operations * cycle_time_ns

# Step 4: Convert the total time from nanoseconds to milliseconds.
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1_000_000

# Step 5: Print the final result, which is the duration in milliseconds.
print(total_time_ms)