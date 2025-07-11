import textwrap

def evaluate_hr4_statements():
    """
    Evaluates multiple-choice statements about the plant gene HR4
    based on scientific literature.
    """

    # Data structure to hold options, their truthfulness, and reasoning.
    # 'score' indicates confidence: 1 for True, 0 for False, 0.9 for technically true but less optimal answer.
    options = {
        'A': {
            'text': "It is an interactor of the actin assembly factor ADF3",
            'score': 0,
            'reasoning': "No scientific evidence supports a direct interaction between HR4 and ADF3. HR4 was identified as interacting with HOP-1, a chaperone-related protein."
        },
        'B': {
            'text': "It contributes to the defense against the broad spectrum of powdery mildew pathogens",
            'score': 1,
            'reasoning': "This is the primary conclusion of the paper that characterized HR4 (Xia et al., 2013, Plant Cell & Environment). The hr4 mutant shows increased susceptibility to multiple species of powdery mildew, confirming its role in broad-spectrum resistance."
        },
        'C': {
            'text': "It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen",
            'score': 0.9,
            'reasoning': "This is also a true statement from the same key paper. HR4-GFP fusion proteins were observed at the Extrahaustorial Membrane (EHM) after infection. However, this describes a mechanistic detail (localization) of its function, while option B describes the overall biological function, making B a more encompassing answer."
        },
        'D': {
            'text': "It regulates the defense modulator PAD4 in plant defense against the Psm.",
            'score': 0,
            'reasoning': "This is incorrect. The hr4 mutant shows no altered resistance to Pseudomonas syringae (Psm), and its salicylic acid pathway, in which PAD4 is a key player, appears largely intact. This suggests HR4's function is specific to powdery mildew and independent of PAD4."
        },
        'E': {
            'text': "HR4 is one of the interactors of the PAD4",
            'score': 0,
            'reasoning': "This is a more specific version of D and is also unsupported by evidence. No direct interaction between HR4 and PAD4 has been reported."
        }
    }

    print("Evaluating the options for the question 'Which of the following is true about HR4?':\n")

    best_answer = None
    max_score = -1

    for option, details in sorted(options.items()):
        print(f"Option {option}: {details['text']}")
        # Use textwrap to format the reasoning nicely
        wrapped_reasoning = textwrap.fill(f"Reasoning: {details['reasoning']}", width=80)
        print(wrapped_reasoning)
        print("-" * 80)
        if details['score'] > max_score:
            max_score = details['score']
            best_answer = option
    
    print(f"\nConclusion: Based on the analysis, the most accurate and encompassing statement is Option {best_answer}.")


# Run the evaluation function
evaluate_hr4_statements()

# Final answer in the required format
# The prompt "Remember in the final code you still need to output each number in the final equation!" is not applicable here.
# We will just print the final letter choice.
print("\nFinal Answer:")
print("<<<B>>>")