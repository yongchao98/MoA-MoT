import random

def play_one_game():
    """
    Simulates one round of the game.

    In each round, we add random numbers (from 0 to 1) to a sum
    until the sum is larger than 1. We then return the count
    of random numbers used.

    Returns:
        int: The number of "throws" needed for the sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1-2: Get a random number between 0 and 1.
        random_value = random.random()
        
        # Step 3: Add to sum and increment the count.
        current_sum += random_value
        num_throws += 1
    return num_throws

def main():
    """
    Runs many simulations and calculates the average number of throws.
    """
    # Step 4: Repeat the process many times to find the average.
    num_simulations = 2000000
    total_throws_recorded = 0

    for _ in range(num_simulations):
        total_throws_recorded += play_one_game()

    # Calculate the average of the recorded numbers.
    average_result = total_throws_recorded / num_simulations

    # Output the final equation with the numbers used.
    print(f"Total throws recorded over all games: {total_throws_recorded}")
    print(f"Total number of games (simulations): {num_simulations}")
    print("---")
    print("The final average is calculated as:")
    print(f"{total_throws_recorded} / {num_simulations} = {average_result}")
    print("\nThis average value converges to the mathematical constant e ≈ 2.71828.")

if __name__ == "__main__":
    main()