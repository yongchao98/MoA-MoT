# Agent C's optimal strategy is to maximize their guaranteed winning probability.
# Through backward induction, we find that C's winning probability V_C,
# as a function of their chosen success probability p_C, has different forms
# over different intervals of p_C in [0, 1].

# A detailed analysis of the subgame perfect equilibrium reveals that for certain
# choices of p_C, agent B will make a move that makes agent A indifferent
# between two strategies. Assuming A breaks the tie against C's favor,
# C's resulting payoff in the region p_C > (sqrt(5)-1)/2 is:
# V_C(p_C) = p_C / (1 + p_C)^2

# Agent C will choose p_C to maximize this function. To find the maximum,
# we can analyze its derivative with respect to p_C.
# Let f(x) = x / (1 + x)^2.
# The derivative f'(x) = (1 - x) / (1 + x)^3.
# Setting the derivative to zero, 1 - x = 0, which gives x = 1.
# This indicates that C's optimal choice for p_C is 1.
p_C = 1.0

# The problem asks for the value of floor(100 * p_C).
# where p_C is the optimal success probability that Agent C should choose.
multiplier = 100
result = int(multiplier * p_C)

# The equation is floor(100 * 1.0)
print(f"The equation to solve is floor({multiplier} * {p_C})")
print(f"The result is {result}")
