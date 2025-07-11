import random

def run_gauss_game():
    """
    Simulates one round of Gauss's game.
    Adds random numbers between 0 and 1 until the sum exceeds 1.
    Returns the count of numbers needed.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number and add it to the sum
        current_sum += random.random()
        num_throws += 1
    # Step 3: The number to be recorded is the count of throws
    return num_throws

def main():
    """
    Runs the simulation multiple times and calculates the average result.
    """
    # Step 4: Reset and repeat the game multiple times
    # We use a relatively small number here to be able to display the full equation.
    num_simulations = 20

    recorded_numbers = []
    for _ in range(num_simulations):
        result = run_gauss_game()
        recorded_numbers.append(result)

    # Calculate the average of the recorded numbers
    total_sum_of_throws = sum(recorded_numbers)
    average = total_sum_of_throws / num_simulations

    # Format the equation string with each recorded number
    equation_str = " + ".join(map(str, recorded_numbers))

    print(f"After {num_simulations} games, the recorded numbers of throws are:")
    print(recorded_numbers)
    print("\nThe average is calculated as follows:")
    # Print the full equation as requested
    print(f"({equation_str}) / {num_simulations} = {average}")

    print("\nThis average is an estimate based on our simulation.")
    print("The theoretical value that the average converges to after a very large")
    print("number of games is Euler's number, 'e'.")

if __name__ == "__main__":
    main()