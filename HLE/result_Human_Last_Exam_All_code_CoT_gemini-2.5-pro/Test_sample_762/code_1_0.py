import math

# Number of loop iterations
iterations = 1000

# Depth of the dependency chain (L -> M -> A -> S)
ops_per_iteration = 4

# Number of parallel execution units
parallel_units = 16

# Calculate the number of batches. Since we can process 'parallel_units' iterations
# at a time for any given stage, the batch size is 16.
# We use math.ceil because even a partial final batch requires a full set of cycles.
num_batches = math.ceil(iterations / parallel_units)

# Each batch takes 'ops_per_iteration' cycles to complete its dependency chain
# because each stage of the chain for the batch occupies all parallel units.
cycles_per_batch = ops_per_iteration

# Calculate the total execution time
shortest_schedule_cycles = num_batches * cycles_per_batch

print(f"Total iterations: {iterations}")
print(f"Operations per iteration (pipeline depth): {ops_per_iteration}")
print(f"Parallel units: {parallel_units}")
print("\n--- Calculation ---")
print(f"Number of batches = ceil(Total Iterations / Parallel Units) = ceil({iterations} / {parallel_units}) = {num_batches}")
print(f"Cycles per batch = Operations per iteration = {cycles_per_batch}")
print(f"Shortest execution schedule (cycles) = Number of batches * Cycles per batch")
print(f"Shortest execution schedule (cycles) = {num_batches} * {cycles_per_batch} = {shortest_schedule_cycles}")