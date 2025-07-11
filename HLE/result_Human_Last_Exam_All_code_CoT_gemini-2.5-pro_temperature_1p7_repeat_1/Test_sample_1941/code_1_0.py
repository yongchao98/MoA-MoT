import sys

def solve_game_of_life_primes():
    """
    Solves the problem by generating special primes and simulating Conway's Game of Life.
    """
    
    # Set this to a lower number (e.g., 50) for a much faster test run.
    # The problem requires 1,000.
    TOTAL_PRIMES_TO_FIND = 1000

    # Step 1: Helper function for primality testing.
    def is_prime(n):
        """Efficiently checks if a number is prime."""
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

    # Step 2: Generate the first 1,000 primes containing only digits 0 and 1.
    print(f"Finding the first {TOTAL_PRIMES_TO_FIND} primes with only 0/1 digits...")
    special_primes = []
    # Start with i=2, whose binary representation '10' gives the integer 10.
    i = 2
    while len(special_primes) < TOTAL_PRIMES_TO_FIND:
        # Generate candidates by interpreting binary strings as decimal integers
        candidate_num = int(bin(i)[2:])
        if is_prime(candidate_num):
            special_primes.append(candidate_num)
        i += 1
        # Provide progress feedback as this step can be slow
        if len(special_primes) % 100 == 0:
             # This flush is for seeing progress in real time
             print(f"Found {len(special_primes)}/{TOTAL_PRIMES_TO_FIND} primes...", file=sys.stdout)
             sys.stdout.flush()


    print(f"Found {len(special_primes)} primes. The first is {special_primes[0]} and the last is {special_primes[-1]}.")
    print("Now, simulating Conway's Game of Life for each pattern...")

    # Step 3: Helper functions for the Game of Life simulation.
    def get_neighbors(cell):
        """Returns the 8 neighbors of a given cell."""
        x, y = cell
        return {
            (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
            (x - 1, y),                 (x + 1, y),
            (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
        }

    def simulate_gol_step(live_cells):
        """Calculates the next generation of live cells."""
        # Find all cells that need to be considered: live cells and their neighbors.
        potential_cells = live_cells.union(*(get_neighbors(c) for c in live_cells))
        
        next_gen_live_cells = set()
        for cell in potential_cells:
            count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
            is_live = cell in live_cells
            
            if is_live and (count == 2 or count == 3):
                next_gen_live_cells.add(cell)
            elif not is_live and count == 3:
                next_gen_live_cells.add(cell)
        return next_gen_live_cells

    def does_pattern_die(number_str, max_generations=500):
        """
        Simulates the pattern for a number string and returns True if it dies out.
        """
        # A pattern is a horizontal line where '1's are live cells.
        live_cells = {(i, 0) for i, digit in enumerate(number_str) if digit == '1'}

        # Keep a history to detect stable patterns or oscillators.
        history = {frozenset(live_cells)}

        for _ in range(max_generations):
            if not live_cells:
                return True  # The pattern has died out.

            live_cells = simulate_gol_step(live_cells)
            
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                return False  # It's a stable or oscillating pattern, so it won't die.
            history.add(frozen_state)

        # Assume it does not die if it runs for max_generations.
        return False

    # Step 4: Count the number of patterns that die out.
    dying_patterns_count = 0
    dying_primes_list = []

    for i, prime_num in enumerate(special_primes):
        if does_pattern_die(str(prime_num)):
            dying_patterns_count += 1
            dying_primes_list.append(prime_num) # In this case we store the number
        
        # Provide progress feedback for the simulation phase
        if (i + 1) % 100 == 0:
            print(f"Simulated {i + 1}/{len(special_primes)} patterns...", file=sys.stdout)
            sys.stdout.flush()

    # Step 5: Output the final answer.
    # The instruction "output each number in the final equation!" is interpreted here
    # as showing the numbers that contribute to the final count.
    print("\nThe prime numbers that form patterns which ultimately die out are:")
    print(dying_primes_list)
    print("\nThe final equation is the count of these numbers.")
    print(f"Total members of C that will ultimately die out = {dying_patterns_count}")
    print(f"\n<<< {dying_patterns_count} >>>")

# Run the solution
solve_game_of_life_primes()