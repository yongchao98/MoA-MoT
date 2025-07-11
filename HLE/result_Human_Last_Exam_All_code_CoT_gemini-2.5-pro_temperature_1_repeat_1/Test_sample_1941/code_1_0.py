import sys

def is_prime(n):
    """
    Checks if a number is prime using an optimized trial division method.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Check only odd divisors up to the square root of n.
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True

def does_pattern_die(number_str):
    """
    Simulates Conway's Game of Life for a 1D pattern represented by the number string.
    Returns True if the pattern eventually has no live cells, False otherwise.
    """
    # Initialize the set of live cells from the input string.
    # The pattern is one-dimensional, so we can fix the y-coordinate at 0.
    live_cells = frozenset({(i, 0) for i, digit in enumerate(number_str) if digit == '1'})

    if not live_cells:
        return True

    history = {live_cells}
    
    # Set a generation limit to catch patterns that might escape to infinity (like gliders)
    # or have extremely long non-repeating lifespans. 500 generations is sufficient for
    # simple starting patterns like these.
    for _ in range(500):
        # Find all cells to consider for the next state: live cells and their neighbors.
        potential_cells = set()
        for x, y in live_cells:
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    potential_cells.add((x + dx, y + dy))

        next_generation = set()
        for x, y in potential_cells:
            # Count the number of live neighbors for each potential cell.
            live_neighbors = 0
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    if dx == 0 and dy == 0:
                        continue
                    if (x + dx, y + dy) in live_cells:
                        live_neighbors += 1

            # Apply the rules of Conway's Game of Life.
            is_currently_alive = (x, y) in live_cells
            if is_currently_alive and live_neighbors in [2, 3]:
                next_generation.add((x, y))
            elif not is_currently_alive and live_neighbors == 3:
                next_generation.add((x, y))
        
        live_cells = frozenset(next_generation)

        if not live_cells:
            return True  # The pattern has died out.

        if live_cells in history:
            return False  # The pattern has entered a stable or oscillating cycle.
        
        history.add(live_cells)

    return False  # Assumed to be stable or ever-growing if it reaches the generation limit.

def solve_game_of_life_primes():
    """
    Finds the first 1,000 prime numbers with only 0/1 digits and counts
    how many of them die out in Conway's Game of Life.
    """
    dying_out_count = 0
    primes_found_count = 0
    
    # We generate number candidates by interpreting integers as binary strings.
    # e.g., i=3 -> '11' -> 11. i=5 -> '101' -> 101.
    # This is an efficient way to generate all numbers containing only digits 0 and 1.
    candidate_generator = 1
    
    while primes_found_count < 1000:
        candidate_generator += 1
        num_str = bin(candidate_generator)[2:]
        
        # Numbers ending in 0 (except 10) are even and not prime. This is a small optimization.
        if len(num_str) > 2 and num_str.endswith('0'):
            continue
            
        num = int(num_str)
        
        if is_prime(num):
            primes_found_count += 1
            if does_pattern_die(num_str):
                dying_out_count += 1
    
    print("From the first 1,000 prime numbers containing only digits 0 and 1:")
    print("The number of patterns that will ultimately die out is:")
    print(dying_out_count)

if __name__ == "__main__":
    solve_game_of_life_primes()
<<<833>>>