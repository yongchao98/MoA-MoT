def find_poet_and_work():
    """
    This function stores and prints information about a specific Spanish poem.
    """
    poet = "Rafael Alberti"
    poem_title = "Ángel de los números"
    collection = "Sobre los ángeles"
    english_poem_title = "Angel of the Numbers"
    english_collection = "Concerning the Angels"
    
    print("The author of these verses is the Spanish poet Rafael Alberti.")
    print("-" * 20)
    print(f"Poet: {poet}")
    print(f"Poem: '{poem_title}' ({english_poem_title})")
    print("\nThe poem is part of a larger work of art, the hugely influential poetry collection:")
    print(f"Work of Art (Collection): '{collection}' ({english_collection})")

# Execute the function to display the answer
find_poet_and_work()