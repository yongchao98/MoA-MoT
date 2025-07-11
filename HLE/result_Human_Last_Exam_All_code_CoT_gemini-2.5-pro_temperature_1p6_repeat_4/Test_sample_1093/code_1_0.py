import random

def simulate_gauss_game():
    """
    Simulates one round of the game.
    Adds random numbers (0 to 1) until the sum exceeds 1.
    Returns the count of numbers added.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Generate a random number between 0 and 1
        current_sum += random.random()
        num_throws += 1
    return num_throws

def main():
    """
    Runs the simulation for a large number of trials and calculates the average.
    """
    # A large number of trials gives a better approximation of the true average.
    num_trials = 2000000
    total_throws = 0

    print(f"Running simulation for {num_trials} trials...")

    # Repeat the game num_trials times and accumulate the results.
    for _ in range(num_trials):
        total_throws += simulate_gauss_game()

    # Calculate the average number of throws.
    average_throws = total_throws / num_trials

    print("\n--- Simulation Results ---")
    print(f"Total number of throws recorded: {total_throws}")
    print(f"Total number of trials: {num_trials}")
    print("\nThe average number of throws converges to the result of the following equation:")
    print(f"{total_throws} / {num_trials} = {average_throws}")
    print("\nThis value is an approximation of Euler's number, e ≈ 2.71828.")

if __name__ == "__main__":
    main()