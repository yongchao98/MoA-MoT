def find_poem_details():
    """
    This function stores and prints the details of a specific Spanish poem.
    """
    poet = "Rafael Alberti"
    artwork = "Pablo Picasso's painting, 'Guernica'"
    poem_verses = """
"No one leaves from here. Nobody. 
Neither the mystic nor the suicidal. 
And it's useless, 
All escape is useless 
(Not even from below
 or from above)."
"""

    print(f"The Spanish poet who wrote these verses is: {poet}")
    print(f"The work of art they were directed to is: {artwork}")
    print("\nThe verses are:")
    print(poem_verses)

if __name__ == "__main__":
    find_poem_details()