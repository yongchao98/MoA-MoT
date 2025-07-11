import random
import math

def simulate_gauss_game():
    """
    Plays one round of the game.
    Returns the number of "throws" (random numbers) needed
    for their sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        random_num = random.random()
        
        # Step 3: Add to sum
        current_sum += random_num
        num_throws += 1
    return num_throws

def main():
    """
    Main function to run the simulation and calculate the average.
    """
    # Step 4: Repeat the process many times
    num_trials = 2000000
    total_throws = 0
    
    print(f"Running simulation for {num_trials} trials...")

    for _ in range(num_trials):
        total_throws += simulate_gauss_game()

    # Calculate the average
    average = total_throws / num_trials

    print("\nSimulation complete.")
    print("The final equation for the average is:")
    # Printing the numbers in the final equation as requested
    print(f"{total_throws} / {num_trials} = {average}")
    
    print(f"\nThe average number of throws converges to: {average}")
    print(f"The exact theoretical answer is e, which is ~{math.e}")

if __name__ == "__main__":
    main()