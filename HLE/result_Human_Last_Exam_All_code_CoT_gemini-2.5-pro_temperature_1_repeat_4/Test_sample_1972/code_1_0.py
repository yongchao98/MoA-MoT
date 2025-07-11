def find_borges_reference():
    """
    This function contains the information about the novel and author
    Jorge Luis Borges described with the specific quote.
    It simulates a lookup and prints the answer.
    """
    
    # A small knowledge base storing the answer
    reference_data = {
        "author": "Juan Carlos Onetti",
        "novel": "El astillero (The Shipyard)",
        "description": "A continuator (and simplifier) of Faulkner",
        "quote": "the intensity of a tiger and the variety that a chess duel can achieve"
    }

    print("Jorge Luis Borges was referring to the following novel and author:")
    print("-" * 60)
    print(f"Author: {reference_data['author']}")
    print(f"Novel: {reference_data['novel']}")
    print("-" * 60)
    print(f"Borges described the author as: '{reference_data['description']}'.")
    print(f"He described the novel's plot as having: '{reference_data['quote']}'.")

# Execute the function to display the information
find_borges_reference()