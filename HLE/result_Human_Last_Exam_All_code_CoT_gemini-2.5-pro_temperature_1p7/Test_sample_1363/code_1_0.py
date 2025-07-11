def find_the_dance():
    """
    This script analyzes a technical ballroom dance question to determine the correct answer.
    It prints the question, the options, and a detailed explanation for the correct choice.
    """
    question = "In which dance is it impossible to overturn a reverse turn without disregarding the technique?"

    options = {
        'A': 'Viennese Waltz',
        'B': 'English Waltz',
        'C': 'European Tango',
        'D': 'Slow Foxtrott',
        'E': 'Quickstep'
    }

    correct_answer = 'C'

    explanation = """
The correct answer is European Tango. Here's the reasoning:

1.  **Swing vs. Staccato:** The Waltz, Viennese Waltz, Foxtrot, and Quickstep are 'swing' dances. They use Contra Body Movement (CBM), sway, and rise & fall to create continuous, flowing rotation. In these dances, a Reverse Turn can be overturned by skilled dancers by increasing the CBM and sway.

2.  **Tango's Unique Technique:** The European Tango is fundamentally different. It is a 'staccato' dance with no sway and no rise & fall. Turns are sharp and executed without the swinging CBM characteristic of the other dances.

3.  **The Impossibility of Overturning:** To add more rotation (overturn) to a Reverse Turn, one would need to introduce sway and a powerful CBM. Doing so in Tango would completely violate its core technical principles (no sway, staccato action). The figure would no longer be recognizable as a Tango turn.

Therefore, it is impossible to overturn a Reverse Turn in the Tango without disregarding the dance's fundamental technique.
"""

    print("--- Question ---")
    print(question)
    print("\n--- Options ---")
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n--- Analysis and Answer ---")
    print(f"The correct option is: {correct_answer}")
    print(explanation)

find_the_dance()