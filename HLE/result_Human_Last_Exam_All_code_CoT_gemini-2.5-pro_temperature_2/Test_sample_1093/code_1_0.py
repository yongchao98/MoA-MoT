import random

def simulate_gauss_game(num_trials):
    """
    Simulates the game described by the user and calculates the average number
    of "needles" needed for the sum of their associated values to exceed 1.
    
    The process is as follows:
    1. A random number between 0 and 1 is generated.
    2. These numbers are added together until their sum is greater than 1.
    3. The count of numbers needed is recorded.
    4. This process is repeated for 'num_trials'.
    5. The final average of all recorded counts is computed.
    """
    
    # A list to store the count from each trial
    recorded_counts = []

    print(f"Running simulation for {num_trials} trials...")

    # Loop for each trial of the game
    for _ in range(num_trials):
        current_sum = 0.0
        num_throws = 0
        # Keep adding random numbers until the sum > 1
        while current_sum <= 1:
            # Generate a random number between 0 and 1
            random_num = random.random()
            current_sum += random_num
            num_throws += 1
        
        # Record the number of throws for this trial
        recorded_counts.append(num_throws)

    # Calculate the average over all trials
    total_throws = sum(recorded_counts)
    average = total_throws / len(recorded_counts)
    
    # The final equation is the sum of all counts divided by the number of trials
    print(f"Sum of all recorded counts: {total_throws}")
    print(f"Number of trials: {len(recorded_counts)}")
    print(f"Average = {total_throws} / {len(recorded_counts)}")
    print(f"\nThe average number of throws converges to: {average}")
    print("This value is an approximation of Euler's number, e (≈ 2.71828).")


# A large number of trials provides a more accurate result.
# Let's run 1,000,000 simulations.
if __name__ == "__main__":
    simulate_gauss_game(1000000)
