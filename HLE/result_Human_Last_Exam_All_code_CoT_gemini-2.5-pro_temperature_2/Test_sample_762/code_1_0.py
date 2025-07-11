import math

def calculate_execution_cycles():
    """
    Calculates the shortest execution schedule for the given C loop on a RISC machine.
    """
    num_iterations = 1000
    ops_per_iteration = 4 # load, multiply, add, store
    # The critical path for a single iteration is LDR -> MUL -> ADD -> STR.
    # With single-cycle latency for each, this path takes 4 cycles.
    latency_per_iteration_chain = 4
    parallel_units = 16

    # We schedule iterations in batches sized to the number of parallel units.
    batch_size = parallel_units
    
    # Calculate the number of full batches.
    num_full_batches = num_iterations // batch_size
    
    # Calculate the number of cycles for these full batches.
    # Each full batch takes 'latency_per_iteration_chain' cycles.
    cycles_for_full_batches = num_full_batches * latency_per_iteration_chain
    
    # Check if there are any remaining iterations that form a final partial batch.
    remaining_iterations = num_iterations % batch_size
    
    cycles_for_last_batch = 0
    if remaining_iterations > 0:
        # The partial batch also takes the full 4-cycle latency to complete.
        cycles_for_last_batch = latency_per_iteration_chain

    # Total cycles is the sum.
    total_cycles = cycles_for_full_batches + cycles_for_last_batch

    print("--- Step-by-step Calculation ---")
    print(f"Number of iterations: {num_iterations}")
    print(f"Operations per iteration: {ops_per_iteration}")
    print(f"Latency of one iteration's dependency chain: {latency_per_iteration_chain} cycles")
    print(f"Number of parallel units: {parallel_units}")
    print("\n--- Scheduling Logic ---")
    print(f"Iterations are processed in batches of {batch_size}.")
    print(f"Number of full batches: {num_iterations} // {batch_size} = {num_full_batches}")
    print(f"Number of remaining iterations in the last batch: {remaining_iterations}")
    print(f"Cycles for full batches: {num_full_batches} batches * {latency_per_iteration_chain} cycles/batch = {cycles_for_full_batches} cycles")
    print(f"Cycles for the final partial batch: {cycles_for_last_batch} cycles")
    print("\n--- Final Equation ---")
    print(f"{cycles_for_full_batches} + {cycles_for_last_batch} = {total_cycles}")

calculate_execution_cycles()