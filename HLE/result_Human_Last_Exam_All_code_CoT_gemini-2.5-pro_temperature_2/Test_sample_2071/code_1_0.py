import sys

def solve_pulp_fiction_dance_query():
    """
    Analyzes the dance scene from 'Pulp Fiction' to identify the specific word
    in the song 'You Never Can Tell' that corresponds to Mia Wallace's
    second exaggerated right leg movement.
    """

    # Scene Information
    movie = "Pulp Fiction"
    song = "You Never Can Tell"
    artist = "Chuck Berry"
    character = "Mia Wallace"

    # Choreography Analysis based on scene observation
    # First Occurrence:
    # During the first chorus, Mia makes the move on the last word of the line.
    first_occurrence_line = "It goes to show you never can tell"
    first_occurrence_word = "tell"
    first_occurrence_leg = "right"

    # Second Occurrence:
    # After a verse, the chorus repeats. Mia makes the same move again.
    second_occurrence_line = "'C'est la vie,' say the old folks, 'it goes to show you never can tell'"
    second_occurrence_word = "tell"
    second_occurrence_leg = "right"

    # The user asks for the word corresponding to the second exaggeration of the right leg.
    answer_word = second_occurrence_word

    print("Task: Identify the word in the song 'You Never Can Tell' from the 'Pulp Fiction' dance scene.")
    print(f"Analyzing the choreography of {character}...\n")
    print(f"First Occurrence of exaggerated right leg movement:")
    print(f"  Line: \"{first_occurrence_line}\"")
    print(f"  Word: \"{first_occurrence_word}\"\n")
    print(f"Second Occurrence of exaggerated right leg movement:")
    print(f"  Line: \"{second_occurrence_line}\"")
    print(f"  Word: \"{second_occurrence_word}\"\n")
    print("------------------------------------------------------------------------------------")
    print(f"The word sung by Chuck Berry when Uma Thurman exaggerates the movement of her right leg for the second time is:")
    print(f"'{answer_word}'")
    print("------------------------------------------------------------------------------------")

if __name__ == '__main__':
    solve_pulp_fiction_dance_query()