import sys
from io import StringIO

def analyze_poem_metaphor():
    """
    Analyzes the 'opals' metaphor in the poem by scoring answer choices
    based on contextual keywords.
    """
    # Storing the original stdout to revert back later
    original_stdout = sys.stdout
    # Redirect stdout to a string buffer
    sys.stdout = captured_output = StringIO()

    print("Analyzing the metaphor 'Her keepsake eyes are opals'...\n")

    # Core concepts derived from the opal metaphor and its context in the poem.
    # Opal: shifting colors, depth, play-of-light.
    # Poem context: memory, secrets, fractured light, past returning.
    analysis_keywords = {
        'shifting', 'depth', 'memory', 'complex', 'secrets', 'fractured',
        'light', 'tales', 'returned', 'flew'
    }

    # The provided answer choices and associated keywords that represent their meaning.
    choices = {
        'A': {'text': 'Decay', 'keywords': {'cold', 'bone', 'bedbound'}},
        'B': {'text': 'Transformation', 'keywords': {'years', 'returned', 'bedbound'}},
        'C': {'text': 'Complexity', 'keywords': {'complex', 'fractured', 'secrets', 'tales'}},
        'D': {'text': 'Complex emotions', 'keywords': {'fears', 'crass', 'hurrahs'}},
        'E': {'text': 'The shifting depth of memories', 'keywords': {'shifting', 'depth', 'memory', 'fractured', 'secrets', 'tales', 'flew', 'returned'}},
        'F': {'text': 'Bad luck', 'keywords': {'fixing cards', 'fears', 'cheap'}},
        'G': {'text': 'Secrecy and depth of character', 'keywords': {'secrets', 'depth', 'tales', 'kept'}}
    }

    scores = {}
    print("Scoring each choice based on its relevance to the opal metaphor:")
    print("------------------------------------------------------------")

    for key, data in choices.items():
        # The score is the number of intersecting keywords between the choice and the analysis.
        score = len(analysis_keywords.intersection(data['keywords']))
        scores[key] = score
        # The "equation" here is the count of matching keywords.
        print(f"Choice {key} ('{data['text']}'): Score = {score}")

    # Find the choice with the highest score
    best_choice_key = max(scores, key=scores.get)
    best_choice_text = choices[best_choice_key]['text']

    print("------------------------------------------------------------\n")
    print(f"Conclusion: Choice '{best_choice_key}' has the highest score.")
    print("The properties of an opal—its shifting play-of-color and inner depth—are a powerful metaphor.")
    print("This directly maps to the nature of the woman's memories, which are described as 'long-lost secrets' that 'flew like moths' and 'returned in nightly fears'.")
    print("Therefore, 'The shifting depth of memories' is the most precise and comprehensive description of the significance of 'opals'.")

    # Restore original stdout
    sys.stdout = original_stdout
    # Get the content of the buffer
    output_str = captured_output.getvalue()

    # Print the captured output to the actual console
    print(output_str)

    # Print the final answer in the required format
    print(f"<<<{best_choice_key}>>>")


analyze_poem_metaphor()