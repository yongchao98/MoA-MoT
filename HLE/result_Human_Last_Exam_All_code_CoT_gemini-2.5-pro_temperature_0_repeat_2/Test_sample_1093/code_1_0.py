import random

def run_gauss_game_trial():
    """
    Simulates one round of the game.

    In each round, we sum random numbers between 0 and 1 until the sum
    is greater than 1.

    Returns:
        int: The number of random numbers (throws) it took for the sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        random_number = random.random()
        current_sum += random_number
        num_throws += 1
    # Step 3: Record the number of throws
    return num_throws

def main():
    """
    Runs a large number of simulations and calculates the average number of throws.
    """
    # Step 4: Reset and repeat many times
    num_trials = 1000000
    total_throws_across_all_trials = 0

    for _ in range(num_trials):
        total_throws_across_all_trials += run_gauss_game_trial()

    # Calculate the final average
    average_result = total_throws_across_all_trials / num_trials

    print(f"Simulation ran for {num_trials} trials.")
    print("The average number of throws needed for the sum to exceed 1 is calculated as:")
    # The final equation as requested
    print(f"{total_throws_across_all_trials} / {num_trials} = {average_result}")
    print("\nAs the number of trials increases, this value converges to Euler's number, e (approx. 2.71828).")

if __name__ == "__main__":
    main()