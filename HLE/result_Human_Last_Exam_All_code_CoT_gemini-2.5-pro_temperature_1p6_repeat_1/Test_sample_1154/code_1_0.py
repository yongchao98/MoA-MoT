# 1. Define the parameters based on the problem statement.
# Test with the highest fault coverage is March RAW, with a complexity of 12N.
complexity_factor = 12
# Number of bits in the RAM.
n_bits = 1000000
# Time for one read/write cycle in nanoseconds.
tc_ns = 5

# 2. Calculate the total time.
# Total time in nanoseconds = complexity * number of bits * time per cycle.
total_time_ns = complexity_factor * n_bits * tc_ns
# Convert nanoseconds to milliseconds (1 ms = 1,000,000 ns).
total_time_ms = total_time_ns / 1_000_000

# 3. Print the equation with the numbers and the final result.
print(f"The chosen test is March RAW with complexity {complexity_factor}N.")
print("The total test time is calculated as: (Complexity Factor * Number of Bits * Time per Cycle) / 1,000,000")
print(f"Test Time (ms) = ({complexity_factor} * {n_bits} * {tc_ns}) / 1000000")
print(f"Result: {total_time_ms} ms")
