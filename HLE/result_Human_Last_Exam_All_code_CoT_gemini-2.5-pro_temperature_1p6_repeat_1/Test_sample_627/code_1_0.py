# A python script to calculate the upper bound for the braid index of the three-twist knot using Vogel's algorithm.

# 1. Define the properties of the standard diagram for the three-twist knot (5_2).
# c is the number of crossings.
c = 5
# s is the number of Seifert circles, found by applying Seifert's algorithm.
s = 2

# 2. State the problem and the formula from Vogel's algorithm.
print("Calculating an upper bound for the braid index of the three-twist knot (5_2) using Vogel's algorithm.")
print(f"The standard diagram of this knot has {c} crossings (c).")
print(f"The number of Seifert circles in this diagram is {s} (s).")
print("\nThe formula for the upper bound derived from Vogel's algorithm is: c - s + 2")

# 3. Calculate and print the result.
upper_bound = c - s + 2
print("\nPlugging in the values:")
print(f"Upper Bound = {c} - {s} + 2")
print(f"Upper Bound = {upper_bound}")
