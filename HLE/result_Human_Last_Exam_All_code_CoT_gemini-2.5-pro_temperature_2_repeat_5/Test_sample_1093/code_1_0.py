import random
import math

def simulate_single_game():
    """
    Simulates one game round.
    Adds random numbers between 0 and 1 until the sum exceeds 1.
    Returns the count of numbers used.
    """
    sum_of_measurements = 0.0
    number_of_throws = 0
    while sum_of_measurements <= 1:
        # Generate a random number between 0 and 1
        sum_of_measurements += random.random()
        number_of_throws += 1
    return number_of_throws

def main():
    """
    Main function to run the simulation and explain the result.
    """
    num_simulations = 2000000 # A large number for a good approximation

    print(f"Running {num_simulations} simulations of the game...")

    total_throws_recorded = 0
    for _ in range(num_simulations):
        total_throws_recorded += simulate_single_game()

    average_result = total_throws_recorded / num_simulations

    print(f"\nThe average number of throws from the simulation is: {average_result}\n")
    print("---------------------------------------------------------")
    print("Theoretical Explanation:")
    print("The average number of throws converges to Euler's number, e.")
    print("This is because the expected number, E[N], is the sum of P(N > k) for k = 0, 1, 2, ...")
    print("P(N > k) is the probability that the sum of k throws is <= 1, which equals 1/k!.")
    print("\nThe final equation is the sum of these probabilities:\n")
    
    terms = [1.0, 1.0, 1/2, 1/6, 1/24, 1/120]
    term_strings = []
    print("P(N > 0) = 1/0! = 1.0")
    print("P(N > 1) = 1/1! = 1.0")
    print("P(N > 2) = 1/2! = 0.5")
    print("P(N > 3) = 1/3! = 0.166666...")
    print("P(N > 4) = 1/4! = 0.041666...")
    print("P(N > 5) = 1/5! = 0.008333...")
    print("... and so on")
    
    print("\nE[N] = 1/0! + 1/1! + 1/2! + 1/3! + ...")
    
    # Building the string for the final equation output
    for i, term_val in enumerate(terms):
        if i < len(terms) -1:
             term_strings.append(f"{term_val:.6f}")
        else:
            term_strings.append("...")

    print(f"Final Equation: {' + '.join(term_strings)} = e ≈ {math.e}")

if __name__ == "__main__":
    main()
