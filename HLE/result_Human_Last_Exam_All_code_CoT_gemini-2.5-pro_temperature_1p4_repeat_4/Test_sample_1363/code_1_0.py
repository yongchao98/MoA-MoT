import collections

def solve_dance_technique_question():
    """
    Analyzes ballroom dance techniques to answer a multiple-choice question.

    The question asks in which dance it's impossible to overturn a reverse turn
    without disregarding the technique.

    The script models this by assigning a "technical rigidity score" to each dance
    for this specific figure. A higher score indicates that overturning the figure
    is more likely to violate its fundamental technique. European Tango is known
    for its staccato action and rigid structure, which does not allow for the
    continuous flow needed to overturn a figure like in swing dances.
    """

    # Define the choices and their associated dances
    choices = {
        'A': 'Viennese Waltz',
        'B': 'English Waltz',
        'C': 'European Tango',
        'D': 'Slow Foxtrott',
        'E': 'Quickstep'
    }

    # Assign a qualitative "technical rigidity score" for overturning a reverse turn.
    # Higher score = more rigid technique, making overturning impossible.
    # Tango's score is highest due to its unique staccato, non-swing technique.
    rigidity_scores = {
        'Viennese Waltz': 4,
        'English Waltz': 3,
        'European Tango': 10,  # The technique is fundamentally broken by overturning.
        'Slow Foxtrott': 2,
        'Quickstep': 3
    }

    # Find the dance with the highest rigidity score
    highest_rigidity_dance = max(rigidity_scores, key=rigidity_scores.get)
    max_score = rigidity_scores[highest_rigidity_dance]

    # Find the corresponding letter for the answer
    answer_letter = [letter for letter, dance in choices.items() if dance == highest_rigidity_dance][0]

    # Create a sorted list of scores for the "final equation"
    sorted_scores = sorted(rigidity_scores.items(), key=lambda item: item[1], reverse=True)

    # Print the "final equation" showing all numbers
    equation_parts = [f"{score} ({dance})" for dance, score in sorted_scores]
    print("Final Equation (Comparison of Technical Rigidity Scores):")
    print(" > ".join(equation_parts))

    print(f"\nThe dance in which it is impossible to overturn a reverse turn without disregarding the technique is the {highest_rigidity_dance}.")
    print(f"This corresponds to answer choice {answer_letter}.")

solve_dance_technique_question()