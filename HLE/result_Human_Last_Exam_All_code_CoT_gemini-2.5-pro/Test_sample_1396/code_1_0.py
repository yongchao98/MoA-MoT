# The problem asks for a single numerical upper bound 'O' for an envy-free
# allocation for 4 agents. While the bound for a 'connected' allocation is
# either infinite (for perfect envy-freeness) or dependent on epsilon,
# a well-known bound exists for a general (not necessarily connected)
# envy-free allocation from the Brams-Taylor protocol.

# The formula for the upper bound on the number of cuts is (n-1)^2,
# where n is the number of agents.

# For n = 4 agents, the final equation is (4 - 1)^2 = 9.
# The numbers that constitute this equation are 4, 1, 2, and 9.

# Define the variables from the equation
num_agents = 4
one = 1
exponent = 2
result = (num_agents - one) ** exponent

# As requested, print each number in the final equation
print(num_agents)
print(one)
print(exponent)
print(result)
