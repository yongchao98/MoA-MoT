import random

def run_simulation():
    """
    This script simulates the game described by J.C.F. Gauss to find the
    average number of random values (between 0 and 1) that must be summed
    to exceed 1.
    """
    num_simulations = 1000000
    total_throws = 0

    # Run the simulation for the specified number of trials
    for _ in range(num_simulations):
        current_sum = 0.0
        throws_in_this_game = 0
        while current_sum <= 1:
            # Step 1 & 2: Get a random number and add it to the sum
            current_sum += random.random()
            throws_in_this_game += 1
        
        # Step 3: Record the number of throws for this game
        total_throws += throws_in_this_game

    # Calculate the average
    average_throws = total_throws / num_simulations

    # Print the final "equation" as requested
    print(f"This simulation was run {num_simulations} times.")
    print("The average number of throws converges to the solution of:")
    print(f"{total_throws} / {num_simulations} = {average_throws}")
    print("\nThis value is an approximation of Euler's number, e (≈2.71828).")

if __name__ == "__main__":
    run_simulation()