def analyze_interest_development():
    """
    Analyzes which phase of Hidi and Renninger's model benefits most
    from concrete, actionable feedback for long-term interest development.
    """

    model_phases = {
        'A': {
            'name': "Triggered Situational Interest",
            'impact_reasoning': "Impact is low. Interest is fleeting and not yet attached to personal value. The student may not be invested enough to act on feedback for long-term development."
        },
        'B': {
            'name': "Maintained Situational Interest",
            'impact_reasoning': "Impact is moderate. Feedback helps engagement in the current supportive context, but since the interest is still externally driven, it may not be enough to foster a self-sustaining individual interest."
        },
        'C': {
            'name': "Emerging Individual Interest",
            'impact_reasoning': "Impact is HIGHLY significant. The student has started to value the topic and engage voluntarily. They are at a critical point where they need to build competence to sustain this new interest. Concrete, next-step feedback directly provides this, boosting their skills and confidence, which is crucial for solidifying a long-term individual interest."
        },
        'D': {
            'name': "Well-Developed Individual Interest",
            'impact_reasoning': "Impact is moderate. The student's interest is already self-sustaining. While they benefit from feedback for skill refinement, it doesn't serve the same foundational interest-building function as it would in an earlier phase."
        },
        'E': {
            'name': "General Statement",
            'impact_reasoning': "This statement is too general. The question asks for the phase with the *most significant long-term impact*, implying a comparative analysis is needed."
        }
    }

    print("Analyzing the impact of concrete, next-step feedback on student interest phases:\n")

    for choice, details in model_phases.items():
        print(f"Choice {choice}: {details['name']}")
        print(f"Analysis: {details['impact_reasoning']}\n")

    print("-" * 50)
    best_choice = 'C'
    print(f"Conclusion: The feedback is most impactful for a student with '{model_phases[best_choice]['name']}'.")
    print("This feedback provides a clear pathway to build competence, which is the key factor in transitioning an emerging interest into a stable, well-developed one.")

# Run the analysis
analyze_interest_development()