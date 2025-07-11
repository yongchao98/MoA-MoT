import math

def permutations(n, k):
    """Calculates the number of permutations P(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // math.factorial(n - k)

def combinations(n, k):
    """Calculates the number of combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# The set of primes between 1 and 5 (for x-moves) are {2, 3}. Wait, the destination is (5,7)
# x-coordinates are {1,2,3,5,7}, y-coordinates are {1,2,3,5,7}
# The intermediate x-primes between 1 and 5 are {2, 3, 7}.
# The intermediate y-primes between 1 and 7 are {2, 3, 5}.
# My plan used this, but the numbers in my plan text above had an error.
# Let's fix. x-primes available between 1 and 5 is {2,3}. 2 primes.
# y-primes available between 1 and 7 is {2,3,5}. 3 primes.
# n_x_primes = 2, namely {2, 3}
# n_y_primes = 3, namely {2, 3, 5}

# Wait, the coordinate 7 is also a prime, so it could be an intermediate step for x, like 1->7->5.
# Let's stick with the interpretation that any prime from the coordinate set can be intermediate.
n_x_choices = 3 # Intermediate x-primes are {2, 3, 7}
n_y_choices = 3 # Intermediate y-primes are {2, 3, 5}

# Case 1: 1 Horizontal, 3 Vertical moves
k1 = 1
moves1_h = permutations(n_x_choices, k1 - 1)
moves1_v = permutations(n_y_choices, 4 - k1 - 1)
skeletons1 = combinations(4, k1)
total1 = skeletons1 * moves1_h * moves1_v

# Case 2: 2 Horizontal, 2 Vertical moves
k2 = 2
moves2_h = permutations(n_x_choices, k2 - 1)
moves2_v = permutations(n_y_choices, 4 - k2 - 1)
skeletons2 = combinations(4, k2)
total2 = skeletons2 * moves2_h * moves2_v

# Case 3: 3 Horizontal, 1 Vertical move
k3 = 3
moves3_h = permutations(n_x_choices, k3 - 1)
moves3_v = permutations(n_y_choices, 4 - k3 - 1)
skeletons3 = combinations(4, k3)
total3 = skeletons3 * moves3_h * moves3_v

total_paths = total1 + total2 + total3

print("Paths with 1 Horizontal, 3 Vertical moves: C(4,1) * P(3,0) * P(3,2) = 4 * 1 * 6 = {}".format(total1))
print("Paths with 2 Horizontal, 2 Vertical moves: C(4,2) * P(3,1) * P(3,1) = 6 * 3 * 3 = {}".format(total2))
print("Paths with 3 Horizontal, 1 Vertical move: C(4,3) * P(3,2) * P(3,0) = 4 * 6 * 1 = {}".format(total3))
print("\nTotal distinct paths = {} + {} + {} = {}".format(total1, total2, total3, total_paths))
