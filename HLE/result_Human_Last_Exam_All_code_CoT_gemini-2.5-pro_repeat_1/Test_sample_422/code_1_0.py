def jtb_classifier(belief_statement, is_true, is_justified):
    """
    Classifies a belief as Knowledge or Belief based on the JTB definition.
    A belief is Knowledge if and only if it is Justified and True.
    """
    print(f"Analyzing Belief: \"{belief_statement}\"")
    print(f"Is it true? {is_true}")
    print(f"Is it justified? {is_justified}")

    if is_justified and is_true:
        epistemic_state = "Knowledge"
    else:
        epistemic_state = "Belief"

    print(f"JTB Conclusion: This is a '{epistemic_state}'.")
    return epistemic_state

def main():
    """
    Demonstrates the two main problems with the JTB definition of knowledge.
    """
    print("--- JTB Knowledge Analysis ---")
    print("Ignoring Gettier problems, let's examine two issues with the JTB model.\n")

    # --- Problem 1: Justified False Beliefs (The Justification-Truth Gap) ---
    print("--- Problem 1: Justified False Beliefs ---")
    print("This problem occurs when a belief is strongly justified but is ultimately false. Justification does not guarantee truth.")
    print("-" * 40)
    belief1 = "The sun revolves around the Earth."
    # For centuries, this was justified by the best available evidence (Ptolemaic system).
    is_true1 = False
    is_justified1 = True
    jtb_classifier(belief1, is_true1, is_justified1)
    print("Explanation: Even a perfectly rational and well-supported belief is not knowledge if it turns out to be false. This highlights a gap between our best justifications and reality.\n")


    # --- Problem 2: The Regress Problem of Justification ---
    print("--- Problem 2: The Regress Problem of Justification ---")
    print("This problem occurs with true beliefs that are not properly justified (e.g., a lucky guess). It raises the question: what is sufficient justification?")
    print("-" * 40)
    belief2 = "My lottery ticket number 12345 will win."
    # Let's say, by pure luck, it does win.
    is_true2 = True
    # The belief was based on a guess, not evidence.
    is_justified2 = False
    jtb_classifier(belief2, is_true2, is_justified2)
    print("Explanation: This is a 'True Belief', but not knowledge. To make it knowledge, it needs justification. But any justification (e.g., 'because a psychic told me') would itself require another justification, leading to a potential infinite regress.")


if __name__ == '__main__':
    main()
