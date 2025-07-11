# Based on the analysis, the number of steps for all three graph classes is Theta(log n).
# f_1(n) for a cycle is Theta(log n).
# f_2(n) for a tree with max degree at most 100 is Theta(log n).
# f_3(n) for a graph with max degree at most 100 is Theta(log n).

# Now we map this complexity to the given categories.
# The function f(n) = Theta(log n) fits the description f(n) = Omega(log n), which is category 9.
# It does not fit into any lower category as log n grows faster than the bounds in categories 1-8.

d1 = 9
d2 = 9
d3 = 9

# The three-digit code is d1d2d3.
print(f"{d1}{d2}{d3}")