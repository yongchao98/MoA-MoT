import random
import math

def simulate_gauss_game():
    """
    Simulates one round of Gauss's game.
    A measurement is a random number between 0 and 1.
    We count how many measurements it takes for their sum to exceed 1.
    """
    current_sum = 0.0
    num_measurements = 0
    while current_sum <= 1:
        # Generate a random number between 0 and 1
        measurement = random.random()
        current_sum += measurement
        num_measurements += 1
    return num_measurements

def main():
    """
    Runs the simulation and prints the explanation and results.
    """
    # --- Explanation ---
    print("The problem is to find the long-term average of the number of measurements")
    print("needed for their sum to exceed 1, where each measurement is a random number")
    print("uniformly distributed between 0 and 1.")
    print("\nThis average, or expected value, converges to Euler's number, e.")
    print("The theoretical value is e ≈ 2.71828.")
    
    print("\nThe value of e is defined by the infinite series:")
    print("e = Σ (from n=0 to ∞) 1/n!")
    
    print("\nHere are the first few terms of the final equation:")
    
    equation_str = "e = "
    value_str = "e = "
    
    # Calculate and print the first 7 terms of the series
    num_terms_to_show = 7
    for n in range(num_terms_to_show):
        term_val = 1 / math.factorial(n)
        equation_str += f"1/{n}!"
        value_str += f"{term_val:.5f}"
        if n < num_terms_to_show - 1:
            equation_str += " + "
            value_str += " + "
    equation_str += " + ..."
    value_str += " + ..."
    
    print(equation_str)
    print(value_str)

    # --- Simulation ---
    num_simulations = 2000000  # A large number for better accuracy
    print(f"\nNow, running a simulation with {num_simulations:,} trials to verify the result.")
    
    total_measurements = 0
    # The sum() function with a generator expression is a concise way to do this
    total_measurements = sum(simulate_gauss_game() for _ in range(num_simulations))
        
    average_measurements = total_measurements / num_simulations
    
    print(f"\nThe average number of measurements from the simulation is: {average_measurements}")

if __name__ == "__main__":
    main()