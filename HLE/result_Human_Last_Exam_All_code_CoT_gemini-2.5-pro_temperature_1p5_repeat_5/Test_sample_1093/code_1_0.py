import random

def run_simulation():
    """
    This function simulates the game described by J.C.F. Gauss to find the
    average number of throws required for the sum to exceed 1.
    The theoretical average converges to Euler's number, e.
    """
    # A large number of simulations will produce a result closer to the theoretical value.
    num_simulations = 2000000
    total_throws = 0

    # Step 4: Repeat the game num_simulations times.
    for _ in range(num_simulations):
        current_sum = 0.0
        throws_in_this_game = 0
        # Step 1-3: Add random numbers until the sum is larger than 1.
        while current_sum <= 1:
            # A random angle divided by 360 gives a uniform random number between 0 and 1.
            random_value = random.random()
            current_sum += random_value
            throws_in_this_game += 1
        
        # Record the number of throws for this game.
        total_throws += throws_in_this_game

    # Calculate the average of the recorded numbers.
    average_result = total_throws / num_simulations

    print("The simulation approximates the average number of throws.")
    print(f"Total simulations run: {num_simulations}")
    print(f"Total throws recorded: {total_throws}")
    print("\nThe average is calculated by dividing total throws by the number of simulations:")
    print(f"{total_throws} / {num_simulations} = {average_result}")
    print("\nAs the number of simulations increases, this average converges to Euler's number, e ≈ 2.71828.")

# Run the simulation and print the results.
run_simulation()
