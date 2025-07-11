import random
import math

def simulate_gauss_game():
    """
    Simulates one round of the game.

    Returns:
        int: The number of throws required for the sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Get a random number between 0 and 1
        throw = random.random()
        current_sum += throw
        num_throws += 1
    return num_throws

def main():
    """
    Runs many simulations and prints the theoretical and empirical results.
    """
    num_simulations = 2000000
    total_throws_recorded = 0

    print("Running simulation...")
    # Step 3 & 4: Repeat the game and record the numbers
    for _ in range(num_simulations):
        total_throws_recorded += simulate_gauss_game()

    # Calculate the average
    empirical_average = total_throws_recorded / num_simulations

    print(f"\n--- Simulation Result ---")
    print(f"After {num_simulations:,} games, the average number of throws was: {empirical_average:.6f}\n")

    # --- Theoretical Explanation ---
    print("--- Theoretical Analysis ---")
    print("The average number of throws converges to a precise mathematical constant.")
    print("The expected value is given by the infinite series:")
    print("\n  Average = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...\n")
    print("Let's look at the numbers in that equation:")

    # Print the first few terms of the series as requested
    equation_str = "  Average = "
    terms_to_show = 6
    term_values = []
    for i in range(terms_to_show):
        term_values.append(f"{1/math.factorial(i):.4f}")
    equation_str += " + ".join(term_values) + " + ..."
    print(equation_str)

    print("\nThis series is the definition of Euler's number, e.")
    print(f"The exact theoretical average is e ≈ {math.e:.6f}")

if __name__ == "__main__":
    main()

<<<2.71828>>>