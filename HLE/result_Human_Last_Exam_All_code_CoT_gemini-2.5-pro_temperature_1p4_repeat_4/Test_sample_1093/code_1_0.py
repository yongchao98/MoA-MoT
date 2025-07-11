import random

def simulate_gauss_game():
    """
    Plays one round of Gauss's game.

    In this game, we repeatedly draw random numbers uniformly distributed
    between 0 and 1 and add them together. The game stops when the
    sum exceeds 1.

    Returns:
        int: The number of random numbers (throws) required to make the sum > 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        throw_value = random.random()
        current_sum += throw_value
        num_throws += 1
    # Step 3: Record the number of throws
    return num_throws

def main():
    """
    Runs a large number of simulations to find the average number of throws.
    """
    # Step 4: Repeat the process many times
    num_simulations = 2000000
    total_throws_recorded = 0

    for _ in range(num_simulations):
        total_throws_recorded += simulate_gauss_game()

    # Calculate the average
    average_result = total_throws_recorded / num_simulations

    print(f"Simulation ran for {num_simulations} games.")
    print("The final average is calculated as follows:")
    
    # The final equation, showing each number in the calculation
    print(f"{total_throws_recorded} / {num_simulations} = {average_result}")
    
    print("\nThis value is an approximation of the mathematical constant 'e'.")

if __name__ == "__main__":
    main()
