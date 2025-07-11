# These are the known values for N(k), the number of sequences of k 
# 90-degree moves that return a Rubik's cube to its solved state.
# N(k) for k=1,2,3 are 0, 12, and 0 respectively.
N4 = 132
N5 = 24
N6 = 1764

# The total number of permutations is calculated using the formula derived from the
# Principle of Inclusion-Exclusion: 132 * N(4) + 12 * N(5) + N(6)
result = 132 * N4 + 12 * N5 + N6

print("The number of desired permutations is given by the expression:")
print("132 * N(4) + 12 * N(5) + N(6)")
print("\nWhere:")
print(f"N(4) = {N4} (sequences of 4 moves returning to solve)")
print(f"N(5) = {N5} (sequences of 5 moves returning to solve)")
print(f"N(6) = {N6} (sequences of 6 moves returning to solve)")
print("\nPlugging in the values, the final equation is:")
print(f"132 * {N4} + 12 * {N5} + {N6} = {result}")
