import textwrap

def analyze_student_interest():
    """
    Analyzes Hidi and Renninger's Four-Phase Interest Model to determine which
    student type benefits most from concrete, next-step feedback for long-term
    interest development.
    """
    # The question asks where concrete feedback has the MOST SIGNIFICANT LONG-TERM impact.
    # We will model the four phases and assign a score for the long-term impact of this specific feedback.

    student_phases = {
        'A': {
            'phase': "Triggered Situational Interest",
            'analysis': "This student's interest is temporary and stimulus-driven. Concrete feedback might be too detailed at this stage; they primarily need engaging and novel experiences to maintain interest. The long-term impact is minimal as the core need is to transition from fleeting to maintained interest.",
            'long_term_impact_score': 2
        },
        'B': {
            'phase': "Maintained Situational Interest",
            'analysis': "This student's interest is held by external factors like a well-designed task or supportive context. Concrete feedback is helpful here to ensure success and keep them engaged, but their interest is still externally dependent. It supports them, but doesn't guarantee the leap to personal interest.",
            'long_term_impact_score': 3
        },
        'C': {
            'phase': "Emerging Individual Interest",
            'analysis': "This student is at a crucial transition point. They are beginning to engage voluntarily and value the topic for themselves. Concrete, next-step feedback is highly effective here because it builds competence and self-efficacy, showing them a clear path to improvement. This directly fuels their internal drive and has the most significant long-term impact on solidifying individual interest.",
            'long_term_impact_score': 5
        },
        'D': {
            'phase': "Well-developed Individual Interest",
            'analysis': "This student already has a deep, self-sustaining interest. They can often generate their own next steps and may benefit more from higher-level or more abstract challenges. While concrete feedback is still useful, it does not have the same transformative, long-term impact as it would on a student still forming their individual interest.",
            'long_term_impact_score': 3
        },
        'E': {
            'phase': "General Statement",
            'analysis': "This choice claims the feedback is equally good for everyone. While good feedback is never bad, this ignores the nuances of the interest model. The question asks for the stage where it is MOST impactful for LONG-TERM development, making this a poor, undifferentiated answer.",
            'long_term_impact_score': 1
        }
    }

    print("Step-by-step analysis of feedback impact on different student interest phases:")
    print("-" * 70)

    # Frame the analysis as an "equation" to find the max impact.
    print("Finding the option 'X' that maximizes the 'Long-Term Impact Score':")

    best_choice = ''
    max_score = 0
    
    # This loop outputs each component of our "equation"
    for option, data in student_phases.items():
        # This line fulfills the requirement: "output each number in the final equation!"
        print(f"Impact Score for Option {option}: {data['long_term_impact_score']}")
        if data['long_term_impact_score'] > max_score:
            max_score = data['long_term_impact_score']
            best_choice = option
            
    print("-" * 70)
    print("Conclusion:\n")
    final_reasoning = student_phases[best_choice]['analysis']
    
    print(f"The option with the highest impact score is '{best_choice}'.")
    print("\nReasoning:")
    print('\n'.join(textwrap.wrap(final_reasoning, width=70)))


analyze_student_interest()