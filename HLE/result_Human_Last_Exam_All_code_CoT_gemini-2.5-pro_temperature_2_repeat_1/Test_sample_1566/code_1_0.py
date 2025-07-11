# The problem is a mathematical one related to topology.
# Through a process of logical deduction based on theorems of continua,
# we conclude that the given properties are contradictory.
#
# 1. The definition of an "end point" implies the continuum must be "chainable".
# 2. Property (1) requires a finite number of end points (>1).
# 3. Property (2) requires exactly two orbits under homeomorphisms. This forces the
#    space to be highly symmetric.
#
# A simple chainable continuum is an arc (a line segment).
# An arc has 2 end points. Let's say the arc is the interval [0,1].
# Number of endpoints N = 2. This satisfies property (1).
#
# For property (2), we check the orbits of the arc [0,1].
# The group of homeomorphisms consists of two functions: identity f(x)=x and reflection f(x)=1-x.
# Orbit of 0: {0, 1}. This is the set of endpoints.
# Orbit of 0.5: {0.5}. This is a second orbit.
# Orbit of 0.25: {0.25, 0.75}. This is a third orbit.
# The arc has infinitely many orbits, so it fails property (2).
#
# More complex chainable continua with more than 2 endpoints are known to have
# endpoints of different topological characters, which prevents them from being
# in the same orbit.
#
# Thus, no such continuum exists. The number of such continua is 0.

# The following code symbolically represents this conclusion.
number_of_continua_found = 0

print(f"Let N be the number of endpoints. We are given 1 < N < infinity.")
print(f"Let O be the number of orbits. We are given O = 2.")
print(f"We tested the arc, which has N = 2. It has O = infinity.")
print(f"We considered chainable continua with N > 2. They fail O = 2 because their endpoints are not all alike.")
print(f"Our step-by-step analysis shows that no continuum satisfies all conditions.")
print(f"The number of such topologically distinct continua is {number_of_continua_found}.")
