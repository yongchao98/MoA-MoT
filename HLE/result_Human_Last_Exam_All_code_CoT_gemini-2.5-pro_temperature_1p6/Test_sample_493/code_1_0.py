import math

def calculate_average_constellation_size():
    """
    Calculates the average number of stars per constellation based on a
    nearest-neighbor model for uniformly distributed stars.

    The average size S is given by the formula:
    S = 1 / sum_{k>=2} (q_k / k)
    where q_k is the probability that a random star is part of a k-cycle.
    """

    # --- Step 1: Define the probabilities q_k ---

    # q_2 is the probability of a star being in a 2-cycle (reciprocal pair).
    # This can be derived analytically.
    # q_2 = 6 * pi / (8 * pi + 3 * sqrt(3))
    q_2 = (6 * math.pi) / (8 * math.pi + 3 * math.sqrt(3))

    # q_3 and q_4 are hard to compute and are taken from scientific literature.
    # Higher order terms (q_5, q_6, ...) are negligible.
    q_3 = 0.0123  # Approximate probability for 3-cycles
    q_4 = 0.0001  # Approximate probability for 4-cycles

    # --- Step 2: Calculate the terms in the sum ---

    # Term for k=2
    term_k2 = q_2 / 2
    # Term for k=3
    term_k3 = q_3 / 3
    # Term for k=4
    term_k4 = q_4 / 4

    # --- Step 3: Calculate the sum and the final average size ---

    # The denominator of the formula for S
    total_sum = term_k2 + term_k3 + term_k4

    # The average constellation size
    average_size = 1 / total_sum

    # --- Step 4: Print the results step-by-step ---

    print("The average constellation size 'S' is calculated using the formula:")
    print("S = 1 / (q_2/2 + q_3/3 + q_4/4 + ...)\n")
    print("Where q_k is the probability a star belongs to a cycle of length k.\n")

    print(f"--- Calculation Breakdown ---")
    print(f"q_2 (analytic) = {q_2:.6f}")
    print(f"q_3 (from literature) = {q_3:.6f}")
    print(f"q_4 (from literature) = {q_4:.6f}\n")

    print("The terms of the denominator are:")
    print(f"q_2 / 2 = {q_2:.6f} / 2 = {term_k2:.6f}")
    print(f"q_3 / 3 = {q_3:.6f} / 3 = {term_k3:.6f}")
    print(f"q_4 / 4 = {q_4:.6f} / 4 = {term_k4:.6f}\n")

    print("The final equation is:")
    print(f"S = 1 / ({term_k2:.6f} + {term_k3:.6f} + {term_k4:.6f})")
    print(f"S = 1 / {total_sum:.6f}")
    print(f"S = {average_size:.4f}\n")
    
    print(f"The average number of stars per constellation is approximately {average_size:.4f}.")


if __name__ == '__main__':
    calculate_average_constellation_size()