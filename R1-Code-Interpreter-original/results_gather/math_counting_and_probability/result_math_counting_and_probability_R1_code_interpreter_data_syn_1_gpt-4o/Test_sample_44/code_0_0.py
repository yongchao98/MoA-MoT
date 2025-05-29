from itertools import combinations_with_replacement

# Number of balls
n = 4
# Number of boxes
k = 2

# Generate all possible distributions of balls into boxes
distributions = list(combinations_with_replacement(range(k), n))

# Convert distributions to counts of balls in each box
unique_distributions = set()
for dist in distributions:
    count = [dist.count(i) for i in range(k)]
    unique_distributions.add(tuple(sorted(count, reverse=True)))

# Since boxes are indistinguishable, we only count unique distributions
print(len(unique_distributions))