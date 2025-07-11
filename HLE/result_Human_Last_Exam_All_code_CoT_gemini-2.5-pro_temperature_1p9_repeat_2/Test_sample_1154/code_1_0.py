# Step 1: Define the parameters for the calculation.

# The selected test is March CL, which has the highest fault coverage among the options.
# Its complexity is 12N, meaning it performs 12 operations for each memory cell.
ops_per_cell = 12

# The size of the RAM is 1Mbit.
# N = 1 Mbit = 1,000,000 bits
num_bits = 1_000_000

# The time for one read or write operation (cycle time).
time_per_op_ns = 5  # in nanoseconds

# Step 2: Calculate the total test time in nanoseconds.
# Total Time = Number of Bits * Operations per Bit * Time per Operation
total_time_ns = num_bits * ops_per_cell * time_per_op_ns

# Step 3: Convert the total time from nanoseconds to milliseconds.
# 1 ms = 1,000,000 ns
ns_in_ms = 1_000_000
total_time_ms = total_time_ns / ns_in_ms

# Step 4: Print the final equation and the result.
# The user requested to see the numbers in the final equation.
print("Chosen test: March CL (12N Complexity)")
print(f"Equation: ({ops_per_cell} ops/bit * {num_bits:,} bits * {time_per_op_ns} ns/op) / {ns_in_ms:,} ns/ms")
print(f"Result: {total_time_ms} ms")
