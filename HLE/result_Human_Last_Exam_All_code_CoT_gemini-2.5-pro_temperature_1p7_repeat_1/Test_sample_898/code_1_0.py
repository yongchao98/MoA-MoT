from sympy import symbols
from sympy.logic.boolalg import Implies, Equivalent, Not, And, Or

def main():
    """
    This program demonstrates that statement C is a representation of the
    distributive law of logic, which is famously not universally true in
    the logic of quantum mechanics.
    """
    # Define the propositions as symbolic variables
    a = symbols('a')
    b = symbols('b')
    c = symbols('c')

    # Statement C: (¬(a ∧ b) → (a ∧ c)) ↔ (a ∧ (b ∨ c))
    # Let's break it down.
    # Left side of the equivalence (↔)
    lhs = Implies(Not(And(a, b)), And(a, c))

    # Right side of the equivalence (↔)
    rhs = And(a, Or(b, c))

    # In classical logic, Implies(P, Q) is equivalent to Or(Not(P), Q)
    # Let's show this transformation for our lhs
    simplified_lhs = Or(Not(Not(And(a, b))), And(a, c))
    # This simplifies further by double negation to Or(And(a,b), And(a,c))
    final_lhs = Or(And(a, b), And(a, c))

    print("Analyzing Statement C: (¬(a ∧ b) → (a ∧ c)) ↔ (a ∧ (b ∨ c))")
    print("-" * 60)
    print("Let's first analyze the sub-expression: (¬(a ∧ b) → (a ∧ c))")
    print("Using the classical logic rule P → Q  ≡  ¬P ∨ Q, this becomes:")
    print("¬(¬(a ∧ b)) ∨ (a ∧ c)")
    print("Applying the law of double negation (¬¬P ≡ P), we get:")
    print("(a ∧ b) ∨ (a ∧ c)")
    print("-" * 60)
    print("Now, substituting this back into statement C, we find that C is equivalent to:")
    print("( (a ∧ b) ∨ (a ∧ c) ) ↔ ( a ∧ (b ∨ c) )")
    print("-" * 60)
    print("This is the distributive law of conjunction over disjunction.")
    print("This law holds in classical logic, but it fails in quantum logic for incompatible")
    print("observables like position and momentum. Its failure is a fundamental, 'observable'")
    print("feature of the quantum world, making it the correct answer.")

if __name__ == "__main__":
    main()
