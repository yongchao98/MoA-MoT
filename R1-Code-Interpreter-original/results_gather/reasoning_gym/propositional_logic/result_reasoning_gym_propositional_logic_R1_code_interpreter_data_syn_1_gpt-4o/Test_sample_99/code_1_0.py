from sympy.logic.boolalg import Implies, And, Or, Not, Equivalent
from sympy.abc import P, Q, R

# Define the premises
premise1 = And(Implies(Q, P), Implies(P, Q))  # (Q → P) ∧ (P → Q)
premise2 = Equivalent(And(R, P), R)           # (R ∧ P) ↔ R
premise3 = Implies(Implies(R, P), Or(R, P))   # (R → P) → (R ∨ P)

# Define the conclusion to verify
conclusion = Implies(R, Q)  # R → Q

# Check if the conclusion is a logical consequence of the premises
is_conclusion_valid = Implies(And(premise1, premise2, premise3), conclusion).simplify()

print(is_conclusion_valid)