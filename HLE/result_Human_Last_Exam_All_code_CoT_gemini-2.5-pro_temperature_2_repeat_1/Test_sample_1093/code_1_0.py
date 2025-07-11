import random

def simulate_one_gauss_game():
    """
    Simulates one round of Gauss's game.

    In each round, we generate random numbers between 0 and 1 and add them
    to a sum. We count how many numbers it takes for the sum to exceed 1.

    Returns:
        int: The count of numbers needed for the sum to be greater than 1.
    """
    current_sum = 0.0
    number_of_throws = 0
    while current_sum <= 1:
        # Step 1-2: Generate a random number between 0 and 1
        random_value = random.random()
        current_sum += random_value
        number_of_throws += 1
    # Step 3: The number is recorded when the sum is larger than 1
    return number_of_throws

def main():
    """
    Main function to run the simulation multiple times and report the average.
    """
    # Step 4: Repeat the game many times to find the average
    num_trials = 2000000
    total_throws = 0

    print(f"Running simulation for {num_trials} trials...")

    for _ in range(num_trials):
        total_throws += simulate_one_gauss_game()

    # Calculate the average
    average_throws = total_throws / num_trials

    print("\nSimulation complete.")
    print("The final average is calculated by dividing the sum of all recorded numbers by the number of games played.")
    print("The final equation is:")
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # The final equation is: total_throws / num_trials = average_throws
    print(f"{total_throws} / {num_trials} = {average_throws}")
    
    print("\nAs the number of trials increases, this average converges to Euler's number (e ≈ 2.71828).")


if __name__ == "__main__":
    main()