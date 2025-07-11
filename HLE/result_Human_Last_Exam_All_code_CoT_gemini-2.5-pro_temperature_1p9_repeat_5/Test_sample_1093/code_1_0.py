import random

def run_one_gauss_game():
    """
    Simulates one round of Gauss's game.

    Returns:
        int: The number of random numbers needed for their sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        random_value = random.random()
        current_sum += random_value
        num_throws += 1
    # Step 3: Record the number
    return num_throws

def main():
    """
    Runs the simulation many times and computes the average number of throws.
    """
    # A large number of trials provides a more accurate average
    num_trials = 2_000_000
    total_throws = 0

    # Step 4: Repeat the game and accumulate the results
    for _ in range(num_trials):
        total_throws += run_one_gauss_game()

    # Calculate the average
    average_throws = total_throws / num_trials

    # Print the final calculation as requested
    print(f"The simulation ran {num_trials:,} times.")
    print(f"The total number of throws recorded was {total_throws:,}.")
    print("\nThe average number of throws is the total throws divided by the number of trials:")
    print(f"{total_throws} / {num_trials} = {average_throws}")
    print("\nAs the number of trials increases, this value converges to Euler's number, e (approximately 2.71828).")

if __name__ == "__main__":
    main()