import textwrap

def present_shakespeare_bible_reference():
    """
    This function presents the researched answer to the user's question.
    It identifies the first book of the Bible, in canonical order,
    mentioned by name in a Shakespeare play.
    """

    book_name = "Job"
    play_name = "The Second Part of King Henry IV"
    act = 1
    scene = 2
    character = "Sir John Falstaff"
    quote = "I am as poor as Job, my lord, but not so patient."

    # Present the final answer clearly
    print(f"The first book of the Bible, in canonical order, to be mentioned by name in a Shakespeare play is the Book of {book_name}.")
    print("\n" + "="*50 + "\n")
    print(f"It is mentioned in the play '{play_name}'.\n")
    print("Here are the specific details of the reference:\n")
    print(f"  Play: {play_name}")
    print(f"  Act: {act}")
    print(f"  Scene: {scene}")
    print(f"  Character: {character}")
    print(f"  Quote: \"{quote}\"")
    print("\n" + "="*50)

# Call the function to display the answer.
present_shakespeare_bible_reference()