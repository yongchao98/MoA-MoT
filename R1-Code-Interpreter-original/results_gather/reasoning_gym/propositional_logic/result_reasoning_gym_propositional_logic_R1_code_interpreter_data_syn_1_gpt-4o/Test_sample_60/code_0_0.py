from sympy.logic.boolalg import Or, And, Not, Implies, Equivalent
from sympy.abc import P, Q, R

# Define the premises
premise1 = Equivalent(Equivalent(R, Q), And(R, R))
premise2 = Or(P, R)
premise3 = Not(And(P, P))
premise4 = Or(Equivalent(P, R), Equivalent(P, R))

# Check the consistency of the premises
consistency_check = And(premise1, premise2, premise3, premise4)

# Print the result
print(consistency_check.simplify())