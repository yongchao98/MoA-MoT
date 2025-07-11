import random

def play_one_round():
    """
    Simulates one round of the game.
    Adds random numbers (0 to 1) until the sum is > 1.
    Returns the count of numbers added.
    """
    sum_of_measurements = 0.0
    number_of_throws = 0
    while sum_of_measurements <= 1:
        # Step 1-2: Get a random number between 0 and 1.
        measurement = random.random()
        sum_of_measurements += measurement
        number_of_throws += 1
    # Step 3: Record the number.
    return number_of_throws

def main():
    """
    Runs the simulation many times to find the average number of throws.
    """
    # Step 4: Repeat the process many times.
    num_simulations = 2000000
    total_throws_recorded = 0

    for _ in range(num_simulations):
        total_throws_recorded += play_one_round()

    # Calculate the average of the recorded numbers.
    average_result = total_throws_recorded / num_simulations

    print(f"The simulation was run {num_simulations} times.")
    print(f"The total number of throws across all rounds was {total_throws_recorded}.")
    print("\nThe average number of throws required for the sum to exceed 1 converges to:")
    # Output the final equation with all its numbers.
    print(f"{total_throws_recorded} / {num_simulations} = {average_result}")
    print("\nAs the number of simulations increases, this average approaches the mathematical constant 'e' (Euler's number).")

if __name__ == '__main__':
    main()
