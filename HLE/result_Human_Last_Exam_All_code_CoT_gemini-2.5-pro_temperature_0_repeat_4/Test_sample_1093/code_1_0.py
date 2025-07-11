import random

def simulate_one_game():
    """
    Simulates one game of adding random numbers until the sum exceeds 1.
    Returns the number of random numbers (throws) required.
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
    Runs a large number of simulations and calculates the average result.
    """
    # A large number of trials provides a better approximation of the true average
    num_trials = 1000000
    total_throws_recorded = 0

    # Run the simulation for the specified number of trials
    for _ in range(num_trials):
        total_throws_recorded += simulate_one_game()

    # Calculate the average
    average = total_throws_recorded / num_trials

    # Print the final calculation as an equation
    print(f"The average is the sum of all recorded numbers divided by the number of games:")
    print(f"{total_throws_recorded} / {num_trials} = {average}")

if __name__ == "__main__":
    main()