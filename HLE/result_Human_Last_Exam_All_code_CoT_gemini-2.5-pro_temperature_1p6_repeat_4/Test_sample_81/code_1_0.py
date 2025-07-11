import textwrap

def analyze_tos_clauses():
    """
    This function analyzes several Terms of Service (TOS) clauses to identify
    which one is likely a contract of adhesion with egregious terms or hides
    a material term, based on the 'Doctrine of Reasonable Expectations'.
    """

    # Explanation of the core legal concept
    print("--- Legal Analysis ---")
    explanation = (
        "All modern website TOS are 'contracts of adhesion' (take-it-or-leave-it). "
        "The key question for enforceability is whether a term is so surprising or unfair "
        "that it violates the 'doctrine of reasonable expectations'. A user is not bound "
        "by a term they couldn't reasonably expect to be in the contract."
    )
    print("\n".join(textwrap.wrap(explanation, width=80)))
    print("-" * 22)

    # Analysis of the options provided
    print("\n--- Clause Analysis ---")
    analysis = {
        'A': "Standard in B2B agreements to prevent competition.",
        'B': "Standard license grant for social media/content platforms.",
        'C': "Late fees are expected for paid services. The rate mentioned, '10.5%' per month, is very high and likely a typo for 1.5% per month, but the concept itself is not surprising.",
        'D': "A standard anti-scraping clause found on most websites.",
        'E': "Contains a highly unusual prohibition in sub-point (vi).",
        'F': "A standard anti-abuse clause, with a modern update regarding AI.",
        'G': "The well-known, though often disliked, 'sponsored content' clause for social media."
    }

    for option, desc in analysis.items():
        print(f"Option {option}: {desc}")

    # Identifying the most problematic clause
    print("\n--- Conclusion ---")
    conclusion = (
        "Option E contains a long list of prohibitions, but one is extremely unusual: "
        "'(vi) research or identify any individuals known to be residing in the State of Illinois, U.S.A.' "
        "A user would not reasonably expect a general service to have a specific prohibition "
        "tied to a single state. This is a material term that is surprising and 'hidden' "
        "by its unusual nature, making it a prime example of a term a court might find unenforceable."
    )
    print("\n".join(textwrap.wrap(conclusion, width=80)))

    # Final Answer
    final_answer = "E"
    print(f"\nThe most likely answer is: {final_answer}")


if __name__ == '__main__':
    analyze_tos_clauses()

<<<E>>>