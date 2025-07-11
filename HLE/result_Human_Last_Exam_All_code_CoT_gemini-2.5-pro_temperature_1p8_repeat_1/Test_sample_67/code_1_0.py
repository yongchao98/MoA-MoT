def explain_modal_logic_domains():
    """
    Explains the validity of the Barcan and Converse Barcan formulas
    in systems with decreasing domains.
    """

    print("--- Analyzing Modal Formulas in Systems with Decreasing Domains ---")
    print("A 'decreasing domain' means if world w' is accessible from w, the set of individuals in w' is a subset of those in w.\n")

    print("Formula 1: The Converse Barcan Formula")
    print("Formula: ∃x □φ(x) → □∃x φ(x)")
    print("Translation:")
    print("  IF there exists an individual 'c' in the current world that has property φ in all accessible worlds,")
    print("  THEN in all accessible worlds, it's true that some individual has property φ.")
    print("\nAnalysis:")
    print("This formula HOLDS in systems with decreasing domains.")
    print("If the premise is true, we have an individual 'c' from the current world w.")
    print("The premise states φ(c) is true in every accessible world w'.")
    print("A standard semantic rule is that for φ(c) to be true in w', 'c' must exist in w'.")
    print("Therefore, in any accessible world w', we have found an individual ('c') that exists and has property φ.")
    print("This directly proves the consequent.")
    print("-" * 20)

    print("\nFormula 2: The Barcan Formula")
    print("Formula: □∃x φ(x) → ∃x □φ(x)")
    print("Translation:")
    print("  IF in every accessible world, there is some individual with property φ,")
    print("  THEN there is a single individual in the current world that has property φ in all accessible worlds.")
    print("\nAnalysis:")
    print("This formula DOES NOT HOLD in systems with decreasing domains.")
    print("A simple counterexample demonstrates this:")
    print("  - Let the current world w0 have domain D(w0) = {'a', 'b'}.")
    print("  - Let an accessible world w1 have domain D(w1) = {'a'}. (This domain is decreasing).")
    print("\n  Let's check the premise '□∃x φ(x)':")
    print("  - At w1, is '∃x φ(x)' true? Yes, if we say φ('a') is true in w1.")
    print("  - At w0, is '∃x φ(x)' true? Yes, if we say φ('b') is true in w0.")
    print("  => So the premise can be made TRUE.\n")

    print("  Now, let's check the consequent '∃x □φ(x)':")
    print("  - Does '□φ(a)' hold? This requires φ('a') to be true in both w0 and w1. We can set φ('a') to be false in w0.")
    print("  - Does '□φ(b)' hold? This requires φ('b') to be true in both w0 and w1. But 'b' does not exist in w1, so φ('b') cannot be true there.")
    print("  => So the consequent is FALSE.\n")

    print("Since the premise can be true while the consequent is false, the formula does not hold.\n")

    print("--- Final Conclusion ---")
    print("The converse Barcan formula holds, but the Barcan formula does not hold in all possible worlds.")


if __name__ == "__main__":
    explain_modal_logic_domains()
