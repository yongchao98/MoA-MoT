# 1. Analysis for f1(n) on a cycle of length n.
# The algorithm's runtime on a cycle is dominated by breaking the cycle into paths and then reducing the length of the longest path.
# The length of the longest path L_k after k steps follows the recurrence L_{k+1} = O(log L_k).
# This leads to a total runtime of f1(n) = Theta(log* n) with high probability.
# According to the provided categories, Theta(log* n) corresponds to category 3.
d1 = 3

# 2. Analysis for f2(n) on any tree on n vertices of degree at most 100.
# The general upper bound for this algorithm on any graph with maximum degree Delta is O(log n).
# For trees with a large constant maximum degree, there are known constructions that provide a lower bound of Omega(log n).
# Therefore, the tight bound is f2(n) = Theta(log n).
# According to the provided categories, Omega(log n) corresponds to category 9.
d2 = 9

# 3. Analysis for f3(n) on any graph on n vertices of degree at most 100.
# This is the general case. As established for trees, the runtime is bounded by O(log n) from above and Omega(log n) from below.
# The tight bound is f3(n) = Theta(log n).
# According to the provided categories, Omega(log n) corresponds to category 9.
d3 = 9

# The final answer is the three-digit number d1d2d3.
final_code = f"{d1}{d2}{d3}"

print(f"The category for f1(n) (cycles) is: {d1}")
print(f"The category for f2(n) (trees) is: {d2}")
print(f"The category for f3(n) (general graphs) is: {d3}")
print(f"The final three-digit code is: {final_code}")