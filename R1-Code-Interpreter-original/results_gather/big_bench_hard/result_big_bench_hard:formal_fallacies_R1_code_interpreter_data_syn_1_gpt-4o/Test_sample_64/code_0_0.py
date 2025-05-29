from sympy.logic.boolalg import Implies
from sympy.abc import A, B, C

# A: Fan of FC Lugano
# B: Devotee of PFC CSKA Moskva
# C: Ex-supporter of SV Zulte Waregem

# Premise: A -> (not B and not C)
premise = Implies(A, ~B & ~C)

# Conclusion: (not B and not C) -> A
conclusion = Implies(~B & ~C, A)

# Check if the premise implies the conclusion
is_valid = Implies(premise, conclusion).simplify()

print(is_valid)