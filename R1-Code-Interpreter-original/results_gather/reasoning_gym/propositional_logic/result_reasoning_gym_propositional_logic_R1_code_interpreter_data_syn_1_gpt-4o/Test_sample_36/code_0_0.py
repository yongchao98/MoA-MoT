from sympy.logic.boolalg import Or, And, Not, Implies, Equivalent
from sympy.abc import Q, S, R, P

# Define the premises
premise1 = Equivalent(And(Q, S), Implies(R, R))
premise2 = Equivalent(R, Q)
premise3 = And(And(S, P), Not(P))

# Check if the premises are consistent
is_consistent = premise1.simplify() and premise2.simplify() and premise3.simplify()

# Print the result
print(is_consistent)