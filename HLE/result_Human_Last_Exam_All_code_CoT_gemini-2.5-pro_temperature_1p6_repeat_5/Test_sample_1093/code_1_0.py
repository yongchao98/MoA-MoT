import random
import math

def simulate_one_game():
    """
    Simulates one round of Gauss's game.
    Returns the number of random numbers summed until the total is > 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Get a random number from U(0, 1)
        current_sum += random.random()
        num_throws += 1
    return num_throws

def main():
    """
    Runs the simulation, explains the theory, and presents the results.
    """
    print("This program calculates the average number of U(0,1) random variables")
    print("that must be summed to be greater than 1.")
    print("-" * 70)

    # Theoretical explanation
    print("THEORETICAL ANALYSIS:")
    print("Let N be the number of throws. The average is the expected value E[N].")
    print("It can be shown that E[N] = sum_{n=0 to infinity} P(N > n).")
    print("The probability P(N > n) is equal to 1/n!.")
    print("Therefore, the average converges to the sum of the series:")
    print("E[N] = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...")
    print("This series is the definition of the mathematical constant e.")
    print("-" * 70)

    # Print each number in the final equation (the series for e)
    print("THEORETICAL CALCULATION (first 10 terms of the series for e):")
    
    theoretical_sum = 0
    # The requirement is to output each number in the final equation
    # We interpret this as showing the terms of the series that sum to e.
    for n in range(10):
        term = 1 / math.factorial(n)
        theoretical_sum += term
        print(f"Term n={n}: 1/{n}! = {term:<10.7f} | Cumulative Sum = {theoretical_sum:.7f}")
    
    print(f"\nThe sum converges to e ≈ {math.e:.7f}")
    print("-" * 70)

    # Simulation
    print("MONTE CARLO SIMULATION:")
    num_simulations = 2000000
    print(f"Running {num_simulations:,} simulations of the game...")
    
    total_throws = 0
    for _ in range(num_simulations):
        total_throws += simulate_one_game()
        
    average_throws = total_throws / num_simulations
    
    print("\n--- SIMULATION RESULT ---")
    print(f"Average number of throws from simulation: {average_throws:.7f}")
    print(f"This is an approximation of the theoretical value e ≈ {math.e:.7f}")
    print("-" * 70)

if __name__ == "__main__":
    main()