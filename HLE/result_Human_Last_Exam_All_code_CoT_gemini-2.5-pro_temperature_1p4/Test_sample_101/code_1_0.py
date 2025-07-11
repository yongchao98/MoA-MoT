import textwrap

def analyze_interest_model():
    """
    Analyzes Hidi and Renninger's Four-Phase Interest Model to determine which
    student type benefits most from concrete, actionable feedback for long-term interest development.
    """

    # Hidi and Renninger's Four-Phase Interest Model
    phases = {
        'A': "Triggered Situational Interest: Interest is temporary and sparked by external stimuli. Concrete feedback is useful, but the impact may not be long-term as the initial interest is not stable.",
        'B': "Maintained Situational Interest: Interest is sustained but relies on supportive external conditions (e.g., a teacher, group). Feedback helps maintain engagement, but the interest is not yet internally driven.",
        'C': "Emerging Individual Interest: The student begins to engage voluntarily and see personal value. This is a critical transition. Concrete feedback provides necessary scaffolding to build competence, which reinforces the student's budding personal interest and helps solidify it for the long term.",
        'D': "Well-Developed Individual Interest: The student is already self-driven and deeply engaged. They are less dependent on feedback for 'immediate next steps' and more capable of directing their own learning.",
        'E': "Universal Effectiveness: This statement claims the feedback is equally good for everyone, which contradicts the premise of a phase-based model where interventions have different effects at different stages."
    }

    print("Analyzing the options based on the Four-Phase Interest Model:")
    print("-" * 60)

    for key, value in phases.items():
        print(f"Option {key}:")
        # Use textwrap to format the explanation nicely
        wrapped_text = textwrap.fill(value, width=60)
        print(wrapped_text)
        print("-" * 60)

    print("Conclusion:")
    print("A student with an 'Emerging Individual Interest' (C) is at the crucial point where they are developing self-driven motivation. Concrete, actionable feedback provides the exact support needed to turn this budding interest into a stable, long-term, individual interest. Therefore, the long-term impact is most significant for this group.")

analyze_interest_model()