def solve_schedule():
    """
    Calculates the shortest execution schedule for the given problem.
    """
    # Number of iterations in the loop
    num_iterations = 1000

    # Number of operations inside a single loop iteration
    # 1. load (t1 = p[i])
    # 2. multiply (t1 * 3)
    # 3. add (... + 2)
    # 4. store (q[i] = t2)
    ops_per_iteration = 4

    # Number of parallel execution units in the RISC machine
    parallel_units = 16

    # Calculate the total number of operations to be executed
    total_operations = num_iterations * ops_per_iteration

    # The shortest possible schedule is limited by the total workload divided by the machine's parallelism.
    # Since the loop iterations are independent, the pipeline can be kept full.
    shortest_schedule_cycles = total_operations / parallel_units

    print(f"Number of iterations: {num_iterations}")
    print(f"Operations per iteration: {ops_per_iteration}")
    print(f"Number of parallel units: {parallel_units}")
    print("-" * 20)
    print(f"Total operations to execute = {num_iterations} * {ops_per_iteration} = {total_operations}")
    print(f"Shortest execution schedule = {total_operations} / {parallel_units} = {int(shortest_schedule_cycles)} cycles")

solve_schedule()