import math

def solve_schedule():
    """
    Calculates the shortest execution schedule for the given C program on a RISC machine.
    """
    num_iterations = 1000
    num_parallel_units = 16
    # The dependency chain for each iteration is Load -> Multiply -> Add -> Store.
    # With single-cycle latency, this chain is 4 cycles long.
    cycles_per_op_chain = 4

    # We can process iterations in blocks, with the block size equal to the number of parallel units.
    # We need to calculate how many such blocks are required to cover all iterations.
    # Since we must complete all iterations, we use ceiling division.
    num_blocks = math.ceil(num_iterations / num_parallel_units)

    # Each block, even a partially filled one, takes 4 cycles to complete due to the dependency chain.
    # The blocks are executed sequentially as each block fully utilizes the machine's resources.
    total_cycles = num_blocks * cycles_per_op_chain

    print(f"Total iterations: {num_iterations}")
    print(f"Parallel execution units: {num_parallel_units}")
    print(f"Dependency chain length (cycles per block): {cycles_per_op_chain}")
    print("")
    print("Calculation:")
    print(f"Number of blocks = ceil(Total Iterations / Parallel Units) = ceil({num_iterations} / {num_parallel_units}) = {num_blocks}")
    print(f"Total cycles = Number of Blocks * Cycles per Block")
    # Final equation output
    print(f"Total cycles = {num_blocks} * {cycles_per_op_chain} = {total_cycles}")

solve_schedule()