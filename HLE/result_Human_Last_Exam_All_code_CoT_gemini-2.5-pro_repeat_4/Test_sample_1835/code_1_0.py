class Predicate:
    """Represents a predicate, like 'is mortal' (F)."""
    def __init__(self, description, symbol):
        self.description = description
        self.symbol = symbol

    def __str__(self):
        return self.symbol

class Quantifier:
    """Represents a quantifier, like 'For all x' (∀x)."""
    def __init__(self, description, symbol):
        self.description = description
        self.symbol = symbol

    def __str__(self):
        return self.symbol

class QuantifiedProposition:
    """Represents a quantified proposition, like ∀x Fx."""
    def __init__(self, quantifier, predicate):
        self.quantifier = quantifier
        self.predicate = predicate
        # The variable 'x' is bound by the quantifier.
        self.bound_variable = 'x'

    def __str__(self):
        # Creates the final logical equation format, e.g., "∀x(F(x))"
        return f"{self.quantifier}({self.predicate}({self.bound_variable}))"

def main():
    print("This script models the reasoning based on the Generality Constraint.\n")

    # Premise 1: You understand 'Fa' ("Socrates is mortal").
    # This implies you possess the general concept of the predicate F.
    predicate_F = Predicate("is mortal", "F")
    print(f"From understanding a proposition like 'Fa', you possess the predicate concept:")
    print(f"  Concept: {predicate_F.description}")
    print(f"  Symbol: {predicate_F}\n")

    # Premise 2: You understand universal quantification.
    quantifier_forall = Quantifier("For all x", "∀x")
    print(f"You also possess the concept of universal quantification:")
    print(f"  Concept: {quantifier_forall.description}")
    print(f"  Symbol: {quantifier_forall}\n")

    # The Generality Constraint implies you can combine possessed concepts.
    # We combine the predicate 'F' with the quantifier '∀x'.
    print("Applying the Generality Constraint, you can combine these concepts to form a new thought.\n")
    final_proposition = QuantifiedProposition(quantifier_forall, predicate_F)

    # Print the resulting proposition in an equation-like format.
    print("Final constructed proposition (your new understanding):")
    # The following line outputs the components and the final "equation".
    print(f"{quantifier_forall} + {predicate_F} => {final_proposition}")

if __name__ == "__main__":
    main()
