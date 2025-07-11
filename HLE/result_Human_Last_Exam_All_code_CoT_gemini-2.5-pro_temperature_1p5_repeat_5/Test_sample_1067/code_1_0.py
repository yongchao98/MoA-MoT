def solve_logic_puzzle():
    """
    This function analyzes the provided logic puzzle and prints the correct answer with its justification.
    The puzzle revolves around resolving a logical contradiction: a dog detects an intruder but does not bark.
    """

    # Define the propositions for the correct answer
    p = "P: The dog detects an intruder."
    q = "Q: The dog barked."
    r = "R: The dog was asleep."

    # State the observed facts which create the paradox
    fact1 = "P is true (The dog detected an intruder)."
    fact2 = "Q is false (The dog did not bark)."

    # The proposed logical structure from Choice C
    premise1 = "(P ∧ ¬R)→Q"
    premise2 = "(¬Q ∧ P)"
    conclusion = "R"

    # Explanation of the logical structure
    print("The correct answer is C. Here is the logical breakdown:")
    print("-" * 50)
    print("Propositions:")
    print(f"  {p}")
    print(f"  {q}")
    print(f"  {r}\n")

    print("Observed Facts:")
    print(f"  {fact1}")
    print(f"  {fact2}\n")

    print("Logical Argument Form:")
    print(f"  Premise 1: {premise1}")
    print(f"  Premise 2: {premise2}")
    print(f"  Conclusion: ∴ {conclusion}\n")

    print("Explanation:")
    print(f"The argument translates to: 'If (the dog detects an intruder AND the dog is NOT asleep), then the dog barks.'")
    print(f"We know that 'the dog did NOT bark AND it DID detect an intruder.'")
    print(f"Therefore, the only logical conclusion is that 'the dog was asleep.'")
    print("\nThis resolves the contradiction by introducing a necessary condition (being awake) for the dog to bark, which fits the facts without violating logic.")
    print("-" * 50)
    print("The final equation is: ")
    print(f"[(P ∧ ¬R)→Q] ∧ (¬Q ∧ P),∴R")


solve_logic_puzzle()
<<<C>>>