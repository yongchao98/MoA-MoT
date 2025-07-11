import math

# Number of loop iterations
iterations = 1000

# Operations per iteration
# 1. Load (t1 = p[i])
# 2. Multiply (t1 * 3)
# 3. Add (... + 2)
# 4. Store (q[i] = t2)
ops_per_iteration = 4

# Total number of operations to execute
total_ops = iterations * ops_per_iteration

# Number of parallel execution units on the RISC machine
parallel_units = 16

# --- Method 1: Theoretical Minimum Time (Resource Limit) ---
# This calculation gives the absolute lower bound on the number of cycles,
# assuming perfect parallelism without considering data dependencies.
min_cycles_theoretical = math.ceil(total_ops / parallel_units)

# --- Method 2: Phased Execution (A valid, simple schedule) ---
# This method respects dependencies by completing all operations of one type
# before starting the next type.

# Time to complete all 1000 operations of a single type (e.g., all loads)
cycles_per_phase = math.ceil(iterations / parallel_units)

# Total time is the sum of the times for the four phases (Load, Mul, Add, Store)
total_cycles_phased = ops_per_iteration * cycles_per_phase

print(f"Number of iterations: {iterations}")
print(f"Operations per iteration: {ops_per_iteration}")
print(f"Total operations: {total_ops}")
print(f"Parallel units: {parallel_units}")
print("-" * 30)
print(f"Theoretical minimum cycles (resource-bound): ceil({total_ops} / {parallel_units}) = {min_cycles_theoretical}")
print(f"Cycles for a single phase (e.g., all loads): ceil({iterations} / {parallel_units}) = {cycles_per_phase}")
print(f"Final Answer: Total cycles for a 4-phase schedule = {ops_per_iteration} * {cycles_per_phase} = {total_cycles_phased}")
