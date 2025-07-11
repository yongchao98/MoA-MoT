# Step 1: Define the given parameters
num_bits = 1000000
# March CL has a complexity of 12N, meaning 12 operations per bit.
complexity = 12
time_per_cycle_ns = 5  # Time for one read/write operation in nanoseconds

# Step 2: Calculate the total number of operations
total_operations = num_bits * complexity

# Step 3: Calculate the total time in nanoseconds
total_time_ns = total_operations * time_per_cycle_ns

# Step 4: Convert the total time to milliseconds
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1000000

# Step 5: Print the explanation and the result
print("Chosen Test: March CL (highest fault coverage with 12N complexity)")
print(f"RAM size (N): {num_bits:,} bits")
print(f"Time per operation (tc): {time_per_cycle_ns} ns")
print("\nCalculating the total test duration:")
print(f"Total time (ms) = ({num_bits} bits * {complexity} ops/bit * {time_per_cycle_ns} ns/op) / 1,000,000 ns/ms")
print(f"Result: {total_time_ms}")