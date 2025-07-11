def find_spouse_name():
    """
    This function stores and retrieves information about the characters
    in Jacek Dukaj's novel "Perfekcyjna niedoskonałość".
    """
    # A dictionary holding facts about the book.
    book_facts = {
        "title": "Perfekcyjna niedoskonałość",
        "main_character": "Adam Zamoyski",
        "spouse": "Angelina"
    }

    # Retrieve the spouse's name.
    spouse_name = book_facts["spouse"]

    # Print the answer.
    print(f"The character that the main character of '{book_facts['title']}' marries is named: {spouse_name}")

# Execute the function to get the answer.
find_spouse_name()