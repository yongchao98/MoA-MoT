import textwrap

def explain_inconsistency_in_dtt():
    """
    Explains the reasoning behind the inconsistency in dependent type theory
    and identifies the problematic axiom.
    """

    print("### Analyzing the Inconsistency in Dependent Type Theory ###")
    print("-" * 60)

    # Step 1: Explain the context (Girard's Paradox)
    explanation_paradox = """
    The problem describes a situation that leads to a famous logical contradiction known as Girard's Paradox. This paradox arises in powerful type theories if they are not carefully restricted. The core of the problem is self-reference: creating a 'type of all types' which then includes itself, leading to a contradiction.
    """
    print("Step 1: Understanding the Paradox")
    print(textwrap.fill(explanation_paradox, width=80))
    print("\n")

    # Step 2: Explain the role of Universes
    explanation_universes = """
    To prevent such paradoxes, dependent type theories (like Coq or Agda) use a hierarchy of universes: Type₀, Type₁, Type₂, etc. A type belongs to a universe, but a universe itself belongs to a higher universe (e.g., the type of Natural Numbers 'Nat' is in Type₀, but Type₀ itself is in Type₁). This stratification prevents a type from containing itself.
    """
    print("Step 2: The Safeguard - Universe Hierarchy")
    print(textwrap.fill(explanation_universes, width=80))
    print("\n")

    # Step 3: Explain how the paradox is constructed
    explanation_construction = """
    Girard's Paradox can be constructed if you can violate this hierarchy. Specifically, if you can create a type `V` in a high universe (e.g., Type₁) that represents all types in a lower universe (e.g., Type₀), and then find a way to create an equivalent copy of `V` back down in Type₀, the system becomes inconsistent. The 'structural recursion' mentioned in the problem is the tool used to build the paradoxical term once this hierarchy is broken.
    """
    print("Step 3: Breaking the Safeguard")
    print(textwrap.fill(explanation_construction, width=80))
    print("\n")

    # Step 4: Identify the culprit axiom
    explanation_axiom = """
    We need to find the axiom that allows this dangerous violation of the universe hierarchy. Let's analyze the choices:

    - (A, B, D, E): Extensionality and proof properties are generally consistent.
    - (F, H, I): Classical logic axioms lead to non-constructivity, not this kind of paradox.
    - (G): The axiom of choice is also consistent.
    - (C) Propositional Resizing: This is the axiom that explicitly deals with relating types across different universes. A powerful resizing axiom allows one to 'resize' a type from a higher universe into a lower one. This is exactly the tool needed to break the universe hierarchy and enable Girard's Paradox.

    Therefore, Propositional Resizing is the axiom that is inconsistent with the other features of the system.
    """
    print("Step 4: Identifying the Problematic Axiom")
    print(textwrap.fill(explanation_axiom, width=80))
    print("-" * 60)
    print("\nFinal Conclusion: The axiom that introduces the inconsistency is Propositional Resizing.")
    print("Answer code: C")


if __name__ == "__main__":
    explain_inconsistency_in_dtt()
    print("\n<<<C>>>")