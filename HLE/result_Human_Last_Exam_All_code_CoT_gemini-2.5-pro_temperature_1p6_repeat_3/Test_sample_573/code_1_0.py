import math

# Define the parameters from the problem
N = 200
k = 13

# The number of available integers to form the sequence is k - (-N) + 1.
# This corresponds to the number of objects in the slice category (Z_N)_{/k}.
num_choices = k + N + 1

print(f"Counting n-simplices for N = {N} and k = {k}.")
print("An n-simplex corresponds to a sequence i_0 <= i_1 <= ... <= i_n,")
print(f"where each i_j is an integer from {-N} to {k}.")
print("The number of such sequences is given by the binomial coefficient C(k+N+n+1, n+1).\n")

# Loop through the required values of n (0 to 5)
for n in range(6):
    # Number of elements in the sequence (the n-simplex)
    num_elements = n + 1
    
    # Calculate the number of n-simplices using the derived formula:
    # C(num_choices + n, n + 1) which is C(k+N+n+1, n+1)
    # The math.comb(a, b) function calculates the binomial coefficient "a choose b".
    num_simplices = math.comb(num_choices + n, num_elements)
    
    # Print the result for the current value of n, showing the formula and the numbers involved.
    # The format is: Number of n-simplices = C(k+N+n+1, n+1) = C(simplified_n, simplified_k) = result
    print(f"For n = {n}:")
    print(f"Number of {n}-simplices = C({k} + {N} + {n} + 1, {n} + 1) = C({num_choices + n}, {num_elements}) = {num_simplices}")
    
<<<139828526724>>>