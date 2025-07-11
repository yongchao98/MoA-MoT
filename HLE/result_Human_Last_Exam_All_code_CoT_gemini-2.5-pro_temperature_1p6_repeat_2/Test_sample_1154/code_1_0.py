# Step 1: Define the parameters for the calculation.

# The test with the highest fault coverage from the list is March CL.
# The complexity of the March CL algorithm is 12N, meaning it performs
# 12 operations (reads or writes) for each bit in the memory.
complexity_factor = 12

# The RAM size is 1Mbit, which is 1,000,000 bits.
ram_size_bits = 1000000

# One read or write cycle takes 5 nanoseconds.
cycle_time_ns = 5

# Step 2: Calculate the total test time.

# Total operations = complexity_factor * ram_size_bits
total_operations = complexity_factor * ram_size_bits

# Total time in nanoseconds = total_operations * cycle_time_ns
total_time_ns = total_operations * cycle_time_ns

# Convert nanoseconds to milliseconds (1 ms = 1,000,000 ns).
total_time_ms = total_time_ns / 1000000

# Step 3: Print the result in a clear format.
print("Chosen Test: March CL (highest fault coverage)")
print(f"Complexity: {complexity_factor}N")
print("-" * 30)
print("Calculation of Test Duration:")
print(f"({complexity_factor} operations/bit * {ram_size_bits:,} bits) * {cycle_time_ns} ns/operation = {total_time_ms:.0f} ms")
print("-" * 30)
print("The total test duration is:")
print(f"{total_time_ms:.0f}")
