import random
import math

def run_simulation(num_games):
    """
    Simulates Gauss's needle game for a specified number of games.
    """
    total_throws = 0
    for _ in range(num_games):
        current_sum = 0.0
        num_throws = 0
        while current_sum <= 1:
            # Step 1 & 2: Generate a random number between 0 and 1
            current_sum += random.random()
            num_throws += 1
        # Step 3: Record the number of throws
        total_throws += num_throws
    
    # Calculate the average
    return total_throws / num_games

def explain_and_solve():
    """
    Explains the theoretical solution and verifies it with a simulation.
    """
    print("The problem is to find the expected number of random numbers U(0,1) that must be summed to exceed 1.")
    print("This value converges to Euler's number, e.\n")
    print("The theoretical value is derived from the sum of an infinite series:")
    
    # Display the equation for e
    # Remember in the final code you still need to output each number in the final equation!
    equation_terms = []
    value_terms = []
    num_display_terms = 8 # Number of terms to show in the equation
    
    for i in range(num_display_terms):
        equation_terms.append(f"1/{i}!")
        value_terms.append(f"{1/math.factorial(i):.4f}")
        
    print("e = " + " + ".join(equation_terms) + " + ...")
    print("e = " + " + ".join(value_terms) + " + ...\n")
    
    print(f"The exact value of e is approximately: {math.e}\n")
    
    print("--- Running Simulation ---")
    num_simulations = 2000000
    print(f"Simulating the game {num_simulations} times...")
    
    average_result = run_simulation(num_simulations)
    
    print(f"The average number of throws from the simulation is: {average_result}")
    
    # Calculate the difference to show how close the simulation is to the theory
    difference = abs(average_result - math.e)
    print(f"The difference between the simulation and the true value of e is: {difference:.6f}")

if __name__ == '__main__':
    explain_and_solve()
