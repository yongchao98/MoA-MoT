import random

def simulate_gauss_game(num_simulations: int):
    """
    Simulates Gauss's needle game a specified number of times and calculates
    the average number of throws required for the sum to exceed 1.

    Args:
        num_simulations: The total number of games to play.
    """
    # This will store the sum of the recorded numbers from all games.
    total_throws_recorded = 0

    # Step 4: Reset and repeat the game many times.
    for _ in range(num_simulations):
        current_sum = 0.0
        num_throws = 0
        # Keep throwing until the sum is larger than 1.
        while current_sum <= 1:
            # Step 1-2: Generate a random number between 0 and 1.
            measurement = random.random()
            current_sum += measurement
            num_throws += 1
        
        # Step 3: Record the number of throws for this game.
        total_throws_recorded += num_throws

    # Calculate the average of the recorded numbers.
    average_result = total_throws_recorded / num_simulations

    print("Simulating the game described by J.C.F. Gauss...")
    print(f"Total games played: {num_simulations}")
    print(f"Sum of recorded throws from all games: {total_throws_recorded}")
    print("\nThe average of the recorded numbers is calculated as:")
    # Output the final equation with all its numbers.
    print(f"{total_throws_recorded} / {num_simulations} = {average_result}")
    print("\nThis average converges to the mathematical constant 'e' (Euler's number), which is approximately 2.71828.")


if __name__ == "__main__":
    # We run a large number of simulations to get a good approximation.
    number_of_games = 2000000
    simulate_gauss_game(number_of_games)
