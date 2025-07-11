import random

def simulate_finite_game(num_boxes, num_runs):
    """
    Simulates a finite version of the box guessing game.

    Args:
      num_boxes: The total number of boxes in the game.
      num_runs: The number of times to simulate the game.

    Setup:
    - There are `num_boxes` boxes.
    - An adversary secretly places a number from 0-9 in each box.
    - Alice's guess will be for the first box (box 0).
    - Alice gets to see the numbers in all other boxes (1 to num_boxes-1).
    - Alice's strategy: She assumes the sum of all numbers is `j (mod 10)`,
      where `j` is a number she chooses randomly from 0-9.
    """
    successful_guesses = 0

    for _ in range(num_runs):
        # 1. Adversary sets the box values
        secret_sequence = [random.randint(0, 9) for _ in range(num_boxes)]
        true_value_of_box_0 = secret_sequence[0]
        true_sum_mod_10 = sum(secret_sequence) % 10

        # 2. Alice observes all boxes except box 0
        sum_of_observed_boxes = sum(secret_sequence[1:])

        # 3. Alice chooses her theory randomly
        # She bets that the total sum modulo 10 is `j`.
        assumed_sum_mod_10 = random.randint(0, 9)

        # 4. Alice forms her guess based on her theory
        # Her guess `g` must satisfy: (g + sum_of_observed_boxes) % 10 = assumed_sum_mod_10
        guess_for_box_0 = (assumed_sum_mod_10 - sum_of_observed_boxes % 10 + 10) % 10

        # 5. Check if her guess was correct
        # This happens only if her random theory was correct.
        if guess_for_box_0 == true_value_of_box_0:
            # This is equivalent to checking if assumed_sum_mod_10 == true_sum_mod_10
            successful_guesses += 1
            
    success_probability = successful_guesses / num_runs
    print(f"Simulating a finite game with {num_boxes} boxes for {num_runs} runs.")
    print(f"Alice's theory about the sum was correct {successful_guesses} times.")
    print(f"Estimated success probability: {success_probability:.4f}")
    print("\nThis strategy yields a probability of 1/10, not 9/10.")

simulate_finite_game(num_boxes=20, num_runs=10000)
print("\nThe argument in the explanation shows that even for the infinite case, a guaranteed probability > 0 is not achievable.")
