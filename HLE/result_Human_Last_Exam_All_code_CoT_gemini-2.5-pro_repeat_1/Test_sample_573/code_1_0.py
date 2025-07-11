import math

# Step 1: Define the given parameters.
N = 200
k = 13

# Step 2: Determine the number of objects in the over-category.
# The objects are the integers j such that k <= j <= N.
# Let M be the count of these objects.
M = N - k + 1

# Print an explanation of the calculation.
print(f"The number of objects in the over-category is M = {N} - {k} + 1 = {M}.")
print("The number of n-simplices is the number of non-decreasing sequences of length n+1")
print(f"chosen from these {M} objects.")
print("This is calculated using the multiset coefficient C(M + n, n + 1).")
print("-" * 50)

# Step 3: Loop through n = 0 to 5 and calculate the number of simplices.
for n in range(6):
    # The number of n-simplices is given by the combination with repetition formula:
    # C(M + (n+1) - 1, n+1) which simplifies to C(M + n, n+1).
    num_simplices = math.comb(M + n, n + 1)
    
    # Print the result for each value of n, showing the formula used.
    print(f"For n = {n}, the number of simplices is C({M} + {n}, {n} + 1) = C({M + n}, {n + 1}) = {num_simplices}")
