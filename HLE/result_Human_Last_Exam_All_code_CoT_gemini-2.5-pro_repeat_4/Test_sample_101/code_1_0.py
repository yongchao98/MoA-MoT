import textwrap

def solve_interest_model_question():
    """
    Analyzes Hidi and Renninger's Four-Phase Interest Model to determine
    which student benefits most from concrete, actionable feedback for
    long-term interest development.
    """

    print("Step 1: Understand the core of the question.")
    print("The goal is to find which student stage sees the most *significant long-term impact on interest development* from receiving concrete feedback with immediate next steps.")
    print("-" * 70)

    print("Step 2: Analyze the four phases of interest development.")
    phases = {
        'A. Triggered Situational Interest': 'Interest is fleeting and sparked by external stimuli. The student is not yet invested enough for feedback to have a long-term impact.',
        'B. Maintained Situational Interest': 'Interest is held by external supports (e.g., an engaging activity). Feedback is part of this external support and helps, but the interest is not yet self-driven.',
        'C. Emerging Individual Interest': 'The student begins to engage voluntarily and values the topic. They are actively trying to build competence and internalize their interest. This is a critical transition.',
        'D. Well-Developed Individual Interest': 'Interest is already strong, deep, and self-sustaining. Feedback helps performance but has less impact on *developing* the core interest.'
    }
    for phase, description in phases.items():
        print(f" - {phase}:\n   {textwrap.fill(description, width=65, initial_indent='   ', subsequent_indent='   ')}\n")
    print("-" * 70)

    print("Step 3: Connect the feedback type to the phases.")
    print("Concrete feedback on 'next steps' provides a clear path to building competence.")
    print("This is most impactful for a student who is *actively trying* to build that competence on their own but may not know how.")
    print("\nThis points directly to the student with 'Emerging Individual Interest.' The feedback empowers their voluntary engagement, builds their skills, and helps solidify their budding interest into a long-term, self-driven one.")
    print("-" * 70)

    print("Step 4: Conclude the analysis.")
    final_answer = 'C'
    print(f"Therefore, the student in phase '{list(phases.keys())[2]}' is the most likely to experience a significant long-term impact.")

    # The prompt requires outputting the final equation, but there is no equation.
    # I will interpret this as outputting the final choice clearly.
    print(f"\nFinal Answer Choice: {final_answer}")
    print(f"<<<{final_answer}>>>")


solve_interest_model_question()