import random

def simulate_one_game():
    """
    Simulates one round of the game.

    Returns:
        int: The number of throws needed for the sum to exceed 1.
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
    Runs the simulation multiple times and calculates the average.
    """
    num_games = 2_000_000  # Number of times to play the game for a good average
    total_throws = 0

    # Run the simulation
    for _ in range(num_games):
        total_throws += simulate_one_game()

    # Calculate the average
    average_result = total_throws / num_games

    # Print the final equation as requested
    print(f"The simulation was run {num_games:,} times.")
    print(f"The total number of throws recorded was {total_throws:,}.")
    print("\nThe average is calculated as:")
    print(f"{total_throws} / {num_games} = {average_result}")
    print("\nThis value converges to Euler's number, e.")


if __name__ == "__main__":
    main()
