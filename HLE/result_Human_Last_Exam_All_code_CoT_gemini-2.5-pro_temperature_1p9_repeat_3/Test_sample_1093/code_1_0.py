import random
import math

def run_simulation(num_trials):
    """
    Performs the Gauss game simulation for a specified number of trials.
    
    Args:
        num_trials (int): The number of times to run the game.

    Returns:
        float: The average number of throws required.
    """
    total_throws_recorded = 0
    # Step 4: Repeat the game many times
    for _ in range(num_trials):
        current_sum = 0.0
        throws_for_this_round = 0
        # Step 1-3: Sum numbers until it's > 1, and record the count
        while current_sum <= 1:
            # Step 1-2: Get a random number from 0-1 (the "angle measurement")
            current_sum += random.random()
            throws_for_this_round += 1
        # Step 3: Record the number of throws for this game
        total_throws_recorded += throws_for_this_round
    
    # Calculate the average over all games
    return total_throws_recorded / num_trials

def main():
    """
    Calculates and explains the theoretical average for Gauss's game
    and verifies it with a simulation.
    """
    print("The average number of throws needed in this game converges to a famous mathematical constant: Euler's number, e.")
    print("The value of e can be calculated with the following infinite series, which represents the sum of probabilities P(N>=k):")
    
    # As requested, here are the numbers in the final equation.
    # The expected value is E[N] = P(N>=1) + P(N>=2) + P(N>=3) + ...
    # This evaluates to: E[N] = 1 + 1 + 1/2! + 1/3! + ...
    # which is exactly the series for e = 1/0! + 1/1! + 1/2! + ...
    
    print("\ne = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + 1/5! + ...")
    
    equation_str = "e = "
    # Printing the first 7 terms of the series and their numerical values
    terms = []
    for i in range(7):
        val = 1 / math.factorial(i)
        terms.append(f"{val:.4f}")
    equation_str += " + ".join(terms) + " + ..."
    
    print(equation_str)
    print(f"\nThe exact theoretical value is e ≈ {math.e}\n")

    print("--- Running Simulation ---")
    # A large number of trials provides a good approximation of the theoretical average
    num_simulations = 2000000
    simulated_average = run_simulation(num_simulations)
    
    print(f"Number of games played: {num_simulations}")
    print(f"Simulated average number of throws: {simulated_average}")


if __name__ == '__main__':
    main()