import random

def run_gauss_game_simulation():
    """
    This script simulates the game described by J.C.F. Gauss to find the
    average number of random numbers (between 0 and 1) that must be
    summed to exceed 1.
    """

    # A higher number of trials leads to a more accurate result.
    num_trials = 1000000

    # This will store the sum of all counts from all trials.
    total_count_sum = 0

    # --- Simulation Loop ---
    # We repeat the game `num_trials` times.
    for _ in range(num_trials):
        # --- Single Game ---
        # Reset the sum and count for a new game.
        current_sum = 0.0
        count = 0
        
        # Keep adding random numbers until the sum is > 1.
        while current_sum <= 1:
            # Generate a random number uniformly distributed between 0 and 1.
            random_num = random.random()
            current_sum += random_num
            count += 1
        
        # Add the count from this game to our total sum.
        total_count_sum += count

    # The average is the total sum of counts divided by the number of games.
    average = total_count_sum / num_trials

    # --- Print the Result ---
    print(f"The simulation ran the game {num_trials:,} times.")
    print("The goal is to find the average number of random draws (from 0 to 1) required for their sum to exceed 1.")
    print("\nThe average is calculated by dividing the total number of draws by the number of games played.")
    
    print("\n--- Final Equation ---")
    print(f"Total number of draws = {total_count_sum}")
    print(f"Number of games = {num_trials}")
    print(f"Average = {total_count_sum} / {num_trials} = {average}")
    print("----------------------")

    print("\nAs the number of trials increases, this average converges to the exact theoretical value, which is Euler's number, e.")
    print("e ≈ 2.718281828459045")

# Run the simulation and print the results
run_gauss_game_simulation()