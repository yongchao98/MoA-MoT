def find_answer():
    """
    This function provides the answer to the user's trivia question
    about G.K. Chesterton.
    """
    # The person who made the statement to G.K. Chesterton.
    # Historical context identifies this person, whom Chesterton called
    # "a very brilliant and distinguished Jewish Zionist".
    person_name = "Dr. Chaim Weizmann"

    # The source where Chesterton recounts this story.
    source_book = "The New Jerusalem"
    year_of_publication = 1920

    # The location of the conversation.
    location = "the Mount of Olives, Jerusalem"

    # Print the detailed answer.
    print(f"The person who spoke with G.K. Chesterton on {location} was {person_name}.")
    print(f"Chesterton recounts this conversation in his {year_of_publication} book, '{source_book}'.")
    print(f"{person_name}, who would later become the first President of Israel, told him:")
    print("'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.'")

# Execute the function to print the answer.
find_answer()