import random
import math

def solve_gauss_game():
    """
    Calculates and simulates the result of Gauss's needle game.
    
    The average number of throws converges to Euler's number, e.
    This function demonstrates this by:
    1. Calculating e from its series expansion.
    2. Simulating the game many times to find the experimental average.
    """

    # Part 1: Theoretical Calculation of e
    print("This problem is a classic probability puzzle. The average of the recorded numbers")
    print("converges to a famous mathematical constant: Euler's number, e.")
    print("\nEuler's number can be defined by the infinite series:")
    print("e = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...")
    print("\nTo show the calculation, here are the first 8 terms of the series and their sum:")
    
    num_terms = 8
    e_approximation = 0.0
    equation_terms = []
    
    for i in range(num_terms):
        term = 1.0 / math.factorial(i)
        e_approximation += term
        # The line below creates one of the "numbers" in the final equation
        equation_terms.append(f"{term:.5f}")

    # This print statement fulfills the "output each number in the final equation" requirement
    equation_str = " + ".join(equation_terms)
    print(f"{equation_str} + ... = {e_approximation:.5f}...")
    print(f"The actual value of e is approximately: {math.e:.5f}\n")


    # Part 2: Simulation of the game
    def run_one_gauss_game():
        """Simulates one round of Gauss's game and returns the number of throws."""
        current_sum = 0.0
        count = 0
        while current_sum <= 1:
            # Step 1-2: Get a random number from a uniform distribution [0,1)
            current_sum += random.random()
            count += 1
        # Step 3: Record the number
        return count

    # Step 4: Repeat and average
    num_simulations = 2000000
    total_count = 0
    
    print(f"Now, simulating the game {num_simulations:,} times...")
    for _ in range(num_simulations):
        total_count += run_one_gauss_game()
        
    average_count = total_count / num_simulations

    print(f"\nThe average number of throws recorded from the simulation is: {average_count:.5f}")

# Run the solution
solve_gauss_game()