import math

# --- Problem Parameters ---
num_iterations = 1000
ops_per_iteration = 4  # Load, Multiply, Add, Store
parallel_units = 16
cycles_per_op_chain = 4 # The length of the dependency chain (Ld -> Mul -> Add -> St)

# --- Calculation ---

# The number of batches of iterations we can run. Since we can process 'parallel_units'
# iterations at the same time, we group the total iterations into batches of this size.
# We use math.ceil because even a partial batch requires processing.
num_batches = math.ceil(num_iterations / parallel_units)
num_full_batches = num_iterations // parallel_units
last_batch_size = num_iterations % parallel_units

# Each full batch takes 'cycles_per_op_chain' cycles to complete because the operations
# for each iteration in the batch must be executed sequentially.
# For example, for a batch of 16 iterations:
# Cycle 1: 16 Loads
# Cycle 2: 16 Multiplies
# Cycle 3: 16 Adds
# Cycle 4: 16 Stores
time_for_full_batches = num_full_batches * cycles_per_op_chain

# The last batch, even if not full, also takes the full 'cycles_per_op_chain'
# because the dependency chain still exists for the remaining iterations.
if last_batch_size > 0:
    time_for_last_batch = cycles_per_op_chain
else:
    time_for_last_batch = 0 # No last partial batch

# The total time is the time for all the full batches plus the time for the final partial batch.
# This is equivalent to num_batches * cycles_per_op_chain
total_cycles = time_for_full_batches + time_for_last_batch

# --- Output ---
print("Step 1: Determine the number of operations per iteration.")
print(f"The loop has 4 dependent operations: Load, Multiply, Add, Store. This takes {cycles_per_op_chain} cycles per iteration.")
print("\nStep 2: Determine the number of batches.")
print(f"With {num_iterations} iterations and {parallel_units} parallel units, we have {num_full_batches} full batches and 1 final batch of {last_batch_size} iterations.")
print(f"Total number of batches = ceil({num_iterations} / {parallel_units}) = {num_batches}")
print("\nStep 3: Calculate total cycles.")
print(f"Each batch requires {cycles_per_op_chain} cycles due to data dependencies.")
print(f"Total Cycles = (Number of Full Batches * Cycles per Batch) + Cycles for Last Batch")
print(f"Total Cycles = ({num_full_batches} * {cycles_per_op_chain}) + {time_for_last_batch}")
print(f"Total Cycles = {num_batches} * {cycles_per_op_chain} = {total_cycles}")
