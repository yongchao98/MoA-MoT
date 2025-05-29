from sympy.logic.boolalg import Or, And, Not, Implies
from sympy.abc import Q, P

# Define the premises
premise1 = Or(Not(Q), Not(P))
premise2 = Not(Or(Q, P))
premise3 = And(Or(Q, P), Implies(P, Q))
premise4 = Q

# Check the consistency of the premises
consistency_check = And(premise1, premise2, premise3, premise4)

# Evaluate the consistency
is_consistent = consistency_check.simplify()

print(is_consistent)