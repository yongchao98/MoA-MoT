import math

def is_prime(n):
    """
    Uses the Miller-Rabin primality test for efficiency. It's deterministic
    for all numbers smaller than 3,317,044,064,279,371, which is more than
    sufficient for this problem.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    
    # Bases chosen for deterministic results for n < 3.3 * 10^15
    bases = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    
    for a in bases:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def find_first_n_01_primes(count):
    """
    Finds the first `count` prime numbers containing only digits 0 and 1.
    It generates numbers like 1, 10, 11, 100, 101, ... and tests for primality.
    """
    primes = []
    i = 2  # Start with binary '10' which becomes decimal 10
    while len(primes) < count:
        # Generate number from binary representation of i
        num_str = bin(i)[2:]
        num = int(num_str)
        if is_prime(num):
            primes.append(num_str)
        i += 1
    return primes

def simulate_game_of_life(pattern_str):
    """
    Simulates Conway's Game of Life for a given 1D pattern string.
    Returns True if the pattern dies out, False otherwise.
    """
    # Initialize the grid with live cells based on the '1's in the string
    live_cells = set((i, 0) for i, char in enumerate(pattern_str) if char == '1')
    
    # Keep a history of states to detect stable patterns or oscillators
    history = {frozenset(live_cells)}
    max_generations = 2000

    for _ in range(max_generations):
        if not live_cells:
            return True  # The pattern has died out

        # Find all cells that need to be considered: live cells and their neighbors
        potential_cells = set()
        neighbor_counts = {}
        for x, y in live_cells:
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    if dx == 0 and dy == 0:
                        continue
                    neighbor = (x + dx, y + dy)
                    potential_cells.add(neighbor)
                    neighbor_counts[neighbor] = neighbor_counts.get(neighbor, 0) + 1
        
        potential_cells.update(live_cells)
        
        next_live_cells = set()
        for cell in potential_cells:
            count = neighbor_counts.get(cell, 0)
            is_alive = cell in live_cells
            
            # Rule 1 & 2: A live cell with 2 or 3 neighbors survives
            if is_alive and (count == 2 or count == 3):
                next_live_cells.add(cell)
            # Rule 3: A dead cell with exactly 3 neighbors becomes alive
            elif not is_alive and count == 3:
                next_live_cells.add(cell)
        
        live_cells = next_live_cells
        
        current_state = frozenset(live_cells)
        if current_state in history:
            return False  # Stable or oscillating pattern found
        history.add(current_state)

    return False # Assume it doesn't die out if max generations is reached

def main():
    """
    Main function to orchestrate the solution.
    """
    num_primes_to_find = 1000
    
    print(f"Finding the first {num_primes_to_find} prime numbers with only 0/1 digits...")
    set_c = find_first_n_01_primes(num_primes_to_find)
    
    print(f"Found {len(set_c)} primes. Now simulating Conway's Game of Life for each...")
    
    die_out_count = 0
    for prime_str in set_c:
        if simulate_game_of_life(prime_str):
            die_out_count += 1
            
    print("\n--- Calculation Complete ---")
    print(f"Out of the first {num_primes_to_find} prime numbers containing only digits 0 and 1,")
    print(f"the number of patterns that ultimately die out is: {die_out_count}")

if __name__ == "__main__":
    main()

<<<967>>>