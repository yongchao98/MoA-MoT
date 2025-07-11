def solve_arabesque_question():
    """
    Analyzes Vaganova arabesques to answer the user's question.

    The question is: which two types of arabesque have the dancer's arm
    extended in front to be on the opposite side as the lifted leg?

    This is equivalent to the forward arm being on the SAME side as the SUPPORTING leg.
    """

    analysis = {
        "First Arabesque": "The arm on the same side as the supporting leg is extended forward. This MATCHES the condition.",
        "Second Arabesque": "The arm on the same side as the lifted leg is extended forward. This DOES NOT match the condition.",
        "Third Arabesque": "Both arms are extended forward. One of these arms is on the same side as the supporting leg. This MATCHES the condition.",
        "Fourth Arabesque": "The arm on the same side as the lifted leg is extended forward (in an effacé position). This DOES NOT match the condition."
    }

    print("Analyzing the four Vaganova arabesque positions:")
    print("-" * 50)
    for position, description in analysis.items():
        print(f"{position}: {description}\n")

    print("-" * 50)
    print("Conclusion: The two arabesques that fit the description are the First and Third.")
    print("This corresponds to answer choice C.")

solve_arabesque_question()