import math
from itertools import product
from collections import deque

def is_prime(n):
    """Checks if a number is prime using an efficient method."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_set_c(count=1000):
    """
    Finds the first 'count' prime numbers containing only digits 0 and 1.
    This is done by generating 0/1-digit numbers and testing for primality.
    """
    print(f"Searching for the first {count} prime numbers with only digits 0 and 1...")
    set_c = []
    # Start with length 2, as '1' is not prime, and we check numbers like '11', '101', etc.
    length = 2
    while len(set_c) < count:
        # Generate all binary-looking strings of a certain length, starting with '1'
        for p in product('01', repeat=length - 1):
            num_str = '1' + ''.join(p)
            num = int(num_str)
            if is_prime(num):
                set_c.append(num)
                if len(set_c) == count:
                    break
        length += 1
    print(f"Found all {count} numbers for set C.")
    return set_c

def will_die_out(pattern_str):
    """
    Simulates Conway's Game of Life for a given 1D pattern string.
    Returns True if the pattern dies out, False otherwise.
    """
    live_cells = set()
    for i, digit in enumerate(pattern_str):
        if digit == '1':
            live_cells.add((i, 0))

    history = set()
    MAX_GENERATIONS = 500  # A safe limit for these patterns

    for generation in range(MAX_GENERATIONS):
        if not live_cells:
            return True  # Pattern died out

        # Use frozenset for hashable states to detect loops
        current_state = frozenset(live_cells)
        if current_state in history:
            return False  # Stable or oscillating pattern found
        history.add(current_state)

        # Find all cells to consider for the next generation
        cells_to_check = set()
        for (x, y) in live_cells:
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    cells_to_check.add((x + dx, y + dy))

        next_live_cells = set()
        for (x, y) in cells_to_check:
            # Count live neighbors
            neighbors = 0
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    if dx == 0 and dy == 0:
                        continue
                    if (x + dx, y + dy) in live_cells:
                        neighbors += 1

            # Apply Game of Life rules
            is_currently_alive = (x, y) in live_cells
            if is_currently_alive and neighbors in [2, 3]:
                next_live_cells.add((x, y))
            elif not is_currently_alive and neighbors == 3:
                next_live_cells.add((x, y))

        live_cells = next_live_cells

    return False  # Assumed to not die out after max generations

def main():
    """Main function to solve the problem."""
    # Step 1: Generate the set C
    set_c = find_set_c(1000)

    # Step 2: Simulate and count how many die out
    print("Simulating Game of Life for each of the 1000 patterns...")
    die_out_count = 0
    for num in set_c:
        pattern_str = str(num)
        if will_die_out(pattern_str):
            die_out_count += 1
    
    print("\n--- Calculation Complete ---")
    print(f"The number of members of set C is: {len(set_c)}")
    print(f"The number of members of C that will ultimately die out is: {die_out_count}")


if __name__ == "__main__":
    main()
<<<26>>>