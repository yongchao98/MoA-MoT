import math

def combinations(n, k):
    """Calculates the number of combinations C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def permutations(n, k):
    """Calculates the number of permutations P(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // math.factorial(n - k)

# Define the set of available prime-like indices within the bounding box
qx = [1, 2, 3, 5]
qy = [1, 2, 3, 5, 7]

# The total number of moves in the path
total_moves = 4

# The number of available intermediate coordinates for x and y
# For x, between 1 and 5, the intermediates are {2, 3}
num_intermediate_x = len(qx) - 2
# For y, between 1 and 7, the intermediates are {2, 3, 5}
num_intermediate_y = len(qy) - 2

# Initialize total paths
total_paths = 0
path_calculations = []
path_results = []

# Case 1: 3 Horizontal, 1 Vertical move
h1, v1 = 3, 1
c1 = combinations(total_moves, h1)
px1 = permutations(num_intermediate_x, h1 - 1)
py1 = permutations(num_intermediate_y, v1 - 1)
paths1 = c1 * px1 * py1
total_paths += paths1
path_calculations.append(f"(C({total_moves},{h1}) * P({num_intermediate_x},{h1-1}) * P({num_intermediate_y},{v1-1}))")
path_results.append(f"({c1} * {px1} * {py1})")

# Case 2: 2 Horizontal, 2 Vertical moves
h2, v2 = 2, 2
c2 = combinations(total_moves, h2)
px2 = permutations(num_intermediate_x, h2 - 1)
py2 = permutations(num_intermediate_y, v2 - 1)
paths2 = c2 * px2 * py2
total_paths += paths2
path_calculations.append(f"(C({total_moves},{h2}) * P({num_intermediate_x},{h2-1}) * P({num_intermediate_y},{v2-1}))")
path_results.append(f"({c2} * {px2} * {py2})")

# Case 3: 1 Horizontal, 3 Vertical moves
h3, v3 = 1, 3
c3 = combinations(total_moves, h3)
px3 = permutations(num_intermediate_x, h3 - 1)
py3 = permutations(num_intermediate_y, v3 - 1)
paths3 = c3 * px3 * py3
total_paths += paths3
path_calculations.append(f"(C({total_moves},{h3}) * P({num_intermediate_x},{h3-1}) * P({num_intermediate_y},{v3-1}))")
path_results.append(f"({c3} * {px3} * {py3})")

# Print the final detailed calculation
print("The total number of distinct Prime Paths is calculated by summing the possibilities for each valid combination of horizontal and vertical moves:")
print("")
print(f"Total Paths = {' + '.join(path_calculations)}")
print(f"            = {' + '.join(path_results)}")
print(f"            = {paths1} + {paths2} + {paths3}")
print(f"            = {total_paths}")
