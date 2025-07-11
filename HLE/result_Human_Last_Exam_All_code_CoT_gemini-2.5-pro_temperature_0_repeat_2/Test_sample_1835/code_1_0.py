class Predicate:
    """Represents a property, like 'is an even number'."""
    def __init__(self, symbol, description):
        self.symbol = symbol
        self.description = description

    def __str__(self):
        return f"Predicate '{self.symbol}' ({self.description})"

class Term:
    """Represents a specific object or number, like '4'."""
    def __init__(self, symbol, description):
        self.symbol = symbol
        self.description = description

    def __str__(self):
        return f"Term '{self.symbol}' ({self.description})"

class Quantifier:
    """Represents a quantifier, like 'For all'."""
    def __init__(self, symbol, description):
        self.symbol = symbol
        self.description = description

    def __str__(self):
        return f"Quantifier '{self.symbol}' ({self.description})"

def main():
    # 1. Define the initial concepts for the proposition "Fa".
    predicate_F = Predicate(symbol="F", description="is an even number")
    term_a = Term(symbol="4", description="the number four")

    print("Step 1: Understand the proposition 'Fa'.")
    print(f"This requires understanding its parts:")
    print(f"  - The predicate F: {predicate_F.description}")
    print(f"  - The term a: {term_a.description}")
    print("-" * 20)

    # 2. Form the proposition "Fa".
    print("Step 2: Form the proposition by combining the concepts.")
    # The final equation here is F(a)
    print(f"Final Equation 1: {predicate_F.symbol}({term_a.symbol})")
    print(f"Meaning: '{term_a.description} {predicate_F.description}'.")
    print("According to the Generality Constraint, understanding this means you have isolated the concept of F.")
    print("-" * 20)

    # 3. Introduce the concept of universal quantification.
    print("Step 3: Assume understanding of universal quantification.")
    quantifier_forall = Quantifier(symbol="∀", description="For all")
    print(f"  - The quantifier: {quantifier_forall.description} (symbol: {quantifier_forall.symbol})")
    print("-" * 20)

    # 4. Recombine the predicate F with the quantifier.
    print("Step 4: Recombine the predicate 'F' with the quantifier to form '∀x Fx'.")
    print("The components of the new proposition are:")
    print(f"  - Component 1 (Quantifier): {quantifier_forall.symbol}")
    print(f"  - Component 2 (Predicate): {predicate_F.symbol}")
    print(f"  - Component 3 (Variable): x")
    print("")
    # The final equation here is ∀x F(x)
    print(f"Final Equation 2: {quantifier_forall.symbol}x {predicate_F.symbol}(x)")
    print(f"Meaning: '{quantifier_forall.description} x, x {predicate_F.description}'.")

if __name__ == "__main__":
    main()