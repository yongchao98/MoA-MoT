import math

# Number of iterations in the loop
num_iterations = 1000

# Number of parallel execution units
parallel_units = 16

# Number of dependent operations (stages) in one iteration (LOAD, MUL, ADD, STORE)
ops_per_iteration = 4

# Calculate the number of batches. Since we can process 'parallel_units' iterations
# at a time, we divide the total iterations by the number of units.
# We use math.ceil because even a partial batch requires processing.
num_batches = math.ceil(num_iterations / parallel_units)

# Each batch, regardless of its size (full 16 or partial 8), must pass through
# the 4 stages of the calculation. Each stage takes 1 cycle.
# So, the time to process one batch is equal to the number of dependent operations.
cycles_per_batch = ops_per_iteration

# The total number of cycles is the number of batches multiplied by the
# number of cycles it takes to complete each batch.
total_cycles = num_batches * cycles_per_batch

# Print the final equation
print(f"The shortest execution schedule can be calculated as follows:")
print(f"Number of batches = ceil({num_iterations} / {parallel_units}) = {num_batches}")
print(f"Cycles per batch = {cycles_per_batch}")
print(f"Total cycles = {num_batches} * {cycles_per_batch} = {total_cycles}")
