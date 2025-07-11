# Step 1: Define the parameters for the calculation.

# The test with the highest fault coverage from the list is March RAW.
# The complexity of March RAW is 14N, meaning 14 operations per bit.
complexity = 14

# N is the number of bits in the RAM.
# 1Mbit = 1,000,000 bits
N_bits = 1000000

# tc is the time for a single read/write cycle in nanoseconds.
tc_ns = 5

# Step 2: Calculate the total time for the test.

# Total operations = complexity * number of bits
total_operations = complexity * N_bits

# Total time in nanoseconds = total operations * time per cycle
total_time_ns = total_operations * tc_ns

# Step 3: Convert the total time from nanoseconds to milliseconds.
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1000000

# Step 4: Print the equation and the final result.
print("Chosen Test: March RAW")
print(f"Complexity: {complexity}N")
print(f"RAM Size (N): {N_bits:,} bits")
print(f"Cycle Time (tc): {tc_ns} ns")
print("\nCalculating total duration...")
print(f"Equation: (Complexity * N_bits * tc_ns) / 1,000,000")
print(f"Calculation: ({complexity} * {N_bits} * {tc_ns}) / 1,000,000 = {total_time_ms} ms")
print(f"\nThe duration of the March RAW test is {total_time_ms} milliseconds.")
