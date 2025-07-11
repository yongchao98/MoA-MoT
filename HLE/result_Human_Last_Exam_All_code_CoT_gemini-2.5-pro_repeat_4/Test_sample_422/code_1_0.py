def demonstrate_jtb_problems():
    """
    This script demonstrates two fundamental problems with the
    Justified True Belief (JTB) definition of knowledge,
    ignoring Gettier-style problems.
    """

    # --- Problem 1: The Regress Problem of Justification ---
    print("--- Problem 1: The Regress Problem of Justification ---")
    print("JTB requires a belief to be 'Justified'. But what justifies the justification?")
    print("This can lead to an infinite regress.\n")

    def is_justified(proposition, knowledge_base, chain):
        """
        A recursive function to check for justification, demonstrating the regress.
        """
        if proposition in chain:
            print(f"Circular reasoning detected! '{proposition}' is already in the justification chain.")
            return False, chain

        chain.append(proposition)

        # A foundational belief that requires no further justification.
        if proposition in knowledge_base.get('axioms', []):
            print(f"Justification for '{proposition}' found: It is a foundational axiom.")
            return True, chain

        justifier = knowledge_base.get('justifications', {}).get(proposition)
        if justifier:
            print(f"To justify '{proposition}', we must first justify its reason: '{justifier}'.")
            return is_justified(justifier, knowledge_base, chain)
        else:
            print(f"Justification chain stops. No foundational axiom found for '{proposition}'.")
            return False, chain

    # A knowledge base where the justification chain never reaches a foundation.
    regress_kb = {
        'justifications': {
            'My car is in the parking lot': 'I remember parking it there',
            'I remember parking it there': 'My memory is reliable today',
            'My memory is reliable today': 'I had a good night\'s sleep',
            'I had a good night\'s sleep': 'My new mattress is effective'
            # The chain regresses without end and without hitting a known axiom.
        },
        'axioms': ['I exist'] # A foundational belief not in the chain.
    }

    print("Attempting to justify the belief: 'My car is in the parking lot'")
    is_p1_knowledge, final_chain = is_justified('My car is in the parking lot', regress_kb, [])
    print(f"\nJustification chain: {' -> '.join(final_chain)}")
    print(f"Conclusion: The justification is not grounded. The 'J' in JTB is not satisfied due to the regress problem.\n")
    print("-" * 50)


    # --- Problem 2: The Inaccessibility of Truth ---
    print("\n--- Problem 2: The Inaccessibility of Truth ---")
    print("JTB requires a belief to be 'True'. But an agent has no direct access to objective truth.")
    print("A well-justified belief can still be false.\n")

    # Let's model objective reality, which the agent cannot see.
    objective_reality = {
        'The sun orbits the Earth': False,
        'The Earth is the center of the universe': False,
        'My keys are in my pocket': True
    }

    # The agent's beliefs and their internal justifications.
    agent_beliefs = {
        'The sun orbits the Earth': 'I see it rise in the east and set in the west every day.',
        'My keys are in my pocket': 'I felt them in my pocket just a minute ago.'
    }

    print("Evaluating beliefs from the agent's internal perspective:")
    for belief, justification in agent_beliefs.items():
        print(f"\nAgent's Belief: '{belief}'")
        print(f"Agent's Justification: '{justification}'")
        # From the agent's perspective, they have a Justified Belief (JB).
        # They assume it is also True (T), and therefore Knowledge.
        print("Agent's Conclusion: This feels like knowledge.")

        # Now, we step outside the agent to check against objective reality.
        is_objectively_true = objective_reality.get(belief)
        print("External Observer's Check:")
        if is_objectively_true:
            print(f"Result: The belief is TRUE. JTB is satisfied. This is Knowledge.")
        else:
            print(f"Result: The belief is FALSE. JTB is NOT satisfied. This is not Knowledge.")

    print("\nConclusion: The agent cannot distinguish a justified TRUE belief from a justified FALSE belief.")
    print("The 'T' in JTB is inaccessible from the agent's perspective, making justification fallible.")


if __name__ == '__main__':
    demonstrate_jtb_problems()