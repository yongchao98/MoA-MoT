def find_common_element():
    """
    Analyzes movie trivia to find a common element in the films of
    Fritz Lang and William Friedkin.
    """

    # Data representing the presence of elements in each director's oeuvre.
    # True means the element is a known feature, False means it is not.
    analysis = {
        'A': {
            'element': 'Aboriginal masks',
            'Lang': False,
            'Friedkin': False,
            'reasoning': 'Not a known motif for either director.'
        },
        'B': {
            'element': 'Magic wands',
            'Lang': False,
            'Friedkin': False,
            'reasoning': 'Not a theme in the work of Lang or Friedkin.'
        },
        'C': {
            'element': 'The first ever cyborgs on screen',
            'Lang': True,
            'Friedkin': False,
            'reasoning': "Lang's 'Metropolis' features the first cyborg, but this is not a theme for Friedkin."
        },
        'D': {
            'element': 'Bugs',
            'Lang': True,
            'Friedkin': True,
            'reasoning': ("Friedkin's 'The Exorcist' contains locusts and his film 'Bug' is about insect paranoia. "
                        "Lang's work, which influenced Friedkin, uses themes of societal plague and infestation (e.g., 'Dr. Mabuse').")
        }
    }

    correct_option = None
    print("Analyzing the options:")
    print("-" * 30)

    for option, data in analysis.items():
        print(f"Option {option}: {data['element']}")
        is_common = data['Lang'] and data['Friedkin']
        print(f"  - Appears in both directors' work? {'Yes' if is_common else 'No'}")
        print(f"  - Reasoning: {data['reasoning']}")
        if is_common:
            correct_option = option
    
    print("-" * 30)
    if correct_option:
        print(f"Conclusion: The correct option is {correct_option} because '{analysis[correct_option]['element']}' appear in the works of both directors.")
    else:
        # This would correspond to option E
        print("Conclusion: No element from A-D is common to both directors.")
        correct_option = 'E'

    # Final answer output, as requested.
    # No equation is relevant here, so we will just output the final letter.
    print("\nFinal Answer Selection:")
    print(f"The final answer is option {correct_option}")

find_common_element()