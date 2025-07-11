import random
import math

def solve_gauss_game():
    """
    This script solves the problem of Gauss's game by simulating it and by
    explaining the theoretical convergence to the mathematical constant 'e'.
    """
    
    # --- Part 1: Simulation ---
    # We will simulate the game many times to find the experimental average.
    
    num_simulations = 2000000

    print(f"Playing Gauss's game {num_simulations:,} times to find the average...")

    def run_one_game():
        """
        Simulates one round of the game: sum random numbers until the total is > 1.
        Returns the number of random numbers needed.
        """
        current_sum = 0.0
        num_throws = 0
        while current_sum <= 1:
            # Step 1 & 2: Generate a random number between 0 and 1
            current_sum += random.random()
            num_throws += 1
        # Step 3: Record the number
        return num_throws

    total_throws_recorded = 0
    # Step 4: Repeat the process
    for _ in range(num_simulations):
        total_throws_recorded += run_one_game()

    simulated_average = total_throws_recorded / num_simulations

    print(f"\nThe average number of throws from the simulation is: {simulated_average}\n")


    # --- Part 2: Theoretical Explanation ---
    # The average converges to Euler's number, 'e'.
    # The expected value E[N] = Σ P(N > k) from k=0 to ∞.
    # For this problem, P(N > k) = 1/k!.
    # So, E[N] = Σ 1/k! = 1/0! + 1/1! + 1/2! + ... which is the series for 'e'.
    
    print("-" * 50)
    print("Theoretical Analysis")
    print("-" * 50)
    print("The average of the recorded numbers converges to Euler's number, e.")
    print("The final equation for 'e' is the infinite sum of the reciprocals of factorials:")
    print("e = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...\n")
    print("Here are the numbers for the first few terms in that equation:")

    num_terms_to_show = 10
    current_sum = 0.0
    
    for i in range(num_terms_to_show):
        term_value = 1 / math.factorial(i)
        current_sum += term_value
        print(f"Term {i}: 1/{i}! = {term_value}")

    print(f"\nThe sum of these {num_terms_to_show} terms is: {current_sum}")
    print(f"As we add more terms, this sum converges to the exact value of e ≈ {math.e}")

# Run the full solution
solve_gauss_game()