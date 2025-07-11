import math

# The number of warriors starting the tournament.
num_warriors = 128

# In a single-elimination tournament, the number of rounds required to find a single
# winner from N participants is log2(N).
num_rounds = math.log2(num_warriors)

# Each round of elimination requires two distinct phases because warriors start in different cities:
# 1. Travel Phase: Half of the warriors travel to meet their opponents. This takes 1 day.
# 2. Fight Phase: The paired-up warriors fight. Since all fights can happen in parallel
#    in different cities, this takes 1 day.
# A warrior cannot do both on the same day, so these phases are sequential for each round.
days_per_round = 2

# The total minimum number of days is the number of rounds multiplied by the days per round.
total_days = num_rounds * days_per_round

print("To find the minimum number of days, we calculate the number of rounds and the time per round.")
print(f"Number of warriors: {num_warriors}")
print(f"Number of elimination rounds: log2({num_warriors}) = {int(num_rounds)} rounds")
print(f"Minimum days required per round (1 for travel + 1 for fighting): {days_per_round} days")
print("\nFinal Calculation:")
print(f"{int(num_rounds)} rounds * {days_per_round} days/round = {int(total_days)} days")