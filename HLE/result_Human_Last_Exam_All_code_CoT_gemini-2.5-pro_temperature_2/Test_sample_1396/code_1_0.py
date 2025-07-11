import math

# The number of agents in the cake-cutting problem.
n = 4

# According to the algorithm by Deng, Qi, and Saberi (2012), the upper bound
# on the query complexity for finding a connected ε-envy-free allocation
# is O(n^2 * log(n/ε)).
# For a given number of agents 'n', the dominant term determining the order
# of the complexity is n^2.
# We will calculate the value of this dominant term for n=4.

result = n ** 2

# The final output will present the equation used to arrive at the answer.
print(f"The upper bound O for a connected ε-envy-free allocation is dominated by the term n^2.")
print(f"For n = {n} agents, the calculation is:")
print(f"{n} * {n} = {result}")
