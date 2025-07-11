import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# The problem states we have a set S of 5 points.
num_points = 5

# The condition is that no proper subcontinuum contains 3 of these points.
# The general theorem states that the largest number of continua in the
# decomposition is C(k, m-1), where k is the number of points and m is the
# size of the forbidden set in a proper subcontinuum.
# Here, k=5 and m=3. So we need to calculate C(5, 3-1) = C(5, 2).
k_for_combination = 3 - 1

# Calculate the result
result = combinations(num_points, k_for_combination)

# The final equation is C(5, 2) = 10.
# The problem asks to print the numbers in the final equation.
print(f"The largest number n is given by the combination formula C(k, m-1).")
print(f"With k={num_points} points and a forbidden size of m={3}, we calculate C({num_points}, {3}-1).")
print(f"So, the final equation is: C({num_points}, {k_for_combination}) = {result}")
print(f"The largest number is {result}.")
