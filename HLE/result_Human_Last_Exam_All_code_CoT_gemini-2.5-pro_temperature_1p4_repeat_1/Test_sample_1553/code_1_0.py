def find_character_info():
    """
    This script contains information about Jacek Dukaj's novel
    "Perfekcyjna niedoskonalosc" to answer the user's question.
    """
    
    # Data representing the book's key relationships
    book_data = {
        "title": "Perfekcyjna niedoskonałość",
        "author": "Jacek Dukaj",
        "main_character": {
            "name": "Adam Zamoyski",
            "spouse": "Angelina"
        }
    }

    # Retrieve the names from the data structure
    main_character_name = book_data["main_character"]["name"]
    spouse_name = book_data["main_character"]["spouse"]

    # Print the final answer
    print(f"In Jacek Dukaj's novel 'Perfekcyjna niedoskonałość', the main character, {main_character_name}, marries {spouse_name}.")

find_character_info()