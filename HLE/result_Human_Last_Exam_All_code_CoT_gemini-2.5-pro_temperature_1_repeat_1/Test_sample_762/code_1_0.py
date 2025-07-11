import math

# --- Problem Parameters ---
# Number of iterations in the loop
total_iterations = 1000
# Number of parallel execution units on the RISC machine
parallel_units = 16
# Number of dependent operations in a single loop iteration (load, mul, add, store)
# Since each has single-cycle latency, this is also the number of cycles per batch.
cycles_per_batch = 4

# --- Calculation ---
# Step 1: Calculate the number of batches needed.
# We can process 'parallel_units' iterations at a time.
# We use math.ceil to round up to the nearest whole number, as even a partial
# batch requires a full processing sequence.
num_batches = math.ceil(total_iterations / parallel_units)

# Step 2: Calculate the total execution time.
# Each batch takes 'cycles_per_batch' to complete due to the dependency chain.
# The batches are processed sequentially because each stage of a batch uses all 16 units.
total_cycles = num_batches * cycles_per_batch

# --- Output ---
print(f"Number of iterations: {total_iterations}")
print(f"Parallel units: {parallel_units}")
print(f"Cycles per batch (due to dependency chain): {cycles_per_batch}")
print("\nStep 1: Calculate the number of batches")
print(f"Number of batches = ceil({total_iterations} / {parallel_units}) = {num_batches}")
print("\nStep 2: Calculate the total cycles")
print(f"Total Cycles = Number of batches * Cycles per batch")
print(f"Total Cycles = {num_batches} * {cycles_per_batch} = {total_cycles}")

print(f"\nThe shortest execution schedule is {total_cycles} cycles.")