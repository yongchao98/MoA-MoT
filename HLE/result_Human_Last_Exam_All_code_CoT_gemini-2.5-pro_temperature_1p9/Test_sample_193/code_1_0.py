# Based on the GTO analysis, the most profitable action with AA is to bet.
# The profitability of a bet is an increasing function of its size in this scenario.
# Therefore, the optimal strategy is to bet the maximum possible amount.
# Since this action is strictly more profitable than any other, it is a pure
# strategy, meaning we should take this action 100% of the time.

# Action details
action = "BET"
# The maximum bet is the effective stack size.
sizing = 1000
frequency_percent = 100

# Per instructions, round both percentages and sizings to the nearest even number.
# In this case, 1000 and 100 are already even, so no change is needed.
final_sizing = int(round(sizing / 2) * 2)
final_frequency = int(round(frequency_percent / 2) * 2)

# The user requested to output each number in the final action string.
# We will print the action, the sizing, and the frequency percentage.
print(f"{action} {final_sizing} {final_frequency}%")