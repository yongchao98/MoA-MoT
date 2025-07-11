import math

def calculate_probability():
    """
    Calculates the probability that for each individual, there exists a unique
    type of item for which they hold strictly more copies than any other individual.
    """
    # Use math.factorial for large integer arithmetic
    fact = math.factorial

    # Total number of ways to distribute the items
    # S = 25! / (5!)^5
    S_numerator = fact(25)
    S_denominator = fact(5)**5
    S = S_numerator // S_denominator

    # Number of favorable distributions, F
    # We focus on the dominant case where each dominant count is 2.
    # Number of 5x5 (0,1)-matrices with row/col sums of 3 is 2040.
    num_dominant_matrices = 2040

    # For each such matrix, the number of ways to form it.
    # Ways(C) = (5!)^5 / (2!)^5 = (120/2)^5 = 60^5
    ways_per_matrix = (fact(5) // fact(2))**5

    # Number of ways for a single fixed permutation of dominant types
    F_perm = num_dominant_matrices * ways_per_matrix

    # Total favorable ways is F_perm multiplied by 5! for all possible permutations.
    F = fact(5) * F_perm

    # The final probability
    P = F / S

    # For the output equation, we construct the expression from its components
    # P = (5! * 2040 * 60^5) / (25! / 5!^5)
    f_val_str = f"({fact(5)} * {num_dominant_matrices} * ({fact(5)} // {fact(2)})^{5})"
    f_val_num = fact(5) * num_dominant_matrices * (fact(5) // fact(2))**5
    
    s_val_str = f"({fact(25)} // {fact(5)}^{5})"
    s_val_num = fact(25) // (fact(5)**5)
    
    # Print the final equation with all numbers
    print(f"P = F / S")
    print(f"F = 5! * (Number of dominant matrices) * (Ways per matrix)")
    print(f"F = {fact(5)} * {num_dominant_matrices} * (({fact(5)} // {fact(2)})**5)")
    print(f"F = {f_val_num}")
    print(f"S = 25! / (5!)^5")
    print(f"S = {s_val_num}")
    print(f"P = {f_val_num} / {s_val_num}")
    print(f"P ≈ {P}")


calculate_probability()
<<<0.0030635543459146525>>>