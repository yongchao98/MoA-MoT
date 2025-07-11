import random

def run_gauss_game_simulation(num_simulations):
    """
    Simulates the game described by J.C.F. Gauss and calculates the
    average number of throws needed for the sum to exceed 1.

    Args:
        num_simulations (int): The number of times to run the game simulation.

    Returns:
        float: The average number of throws.
    """
    total_throws = 0

    # Loop to run the simulation the specified number of times
    for _ in range(num_simulations):
        
        # This block simulates one game according to the rules
        sum_of_measurements = 0.0
        throws_in_this_game = 0
        while sum_of_measurements <= 1:
            # Step 1-2: Generate a random number between 0 and 1
            random_value = random.random()
            sum_of_measurements += random_value
            throws_in_this_game += 1
        
        # Step 3: Add the recorded number to our total
        total_throws += throws_in_this_game

    # Calculate the average of the recorded numbers over all simulations
    average = total_throws / num_simulations
    
    # The final equation is the one used to calculate the average.
    # The following print statements output each number in that final equation.
    print(f"Total throws across all games: {total_throws}")
    print(f"Number of games simulated: {num_simulations}")
    print(f"The final average is calculated by the equation: {total_throws} / {num_simulations}")
    
    return average

if __name__ == '__main__':
    # Set the number of simulations. A higher number yields a more accurate result.
    simulations_to_run = 2000000
    
    average_result = run_gauss_game_simulation(simulations_to_run)
    
    print("\nThe average of the recorded numbers converges to:")
    print(average_result)
    print("\nThis value is an approximation of the mathematical constant 'e' (approx. 2.71828).")
