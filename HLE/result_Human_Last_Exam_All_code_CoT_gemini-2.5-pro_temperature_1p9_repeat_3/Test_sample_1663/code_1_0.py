def find_poet_and_work():
    """
    Identifies the poet and the work of art associated with the provided verses
    and prints the information.
    """
    poet = "Rafael Alberti"
    work_of_art = "Sobre los ángeles (Concerning the Angels)"
    poem_title = "Ángel de los números (Angel of the Numbers)"
    verses = """No one leaves from here. Nobody.
Neither the mystic nor the suicidal.
And it's useless,
All escape is useless
(Not even from below
or from above)."""

    print(f"The Spanish poet who wrote these verses is {poet}.")
    print(f"They are from the poem '{poem_title}' which is part of his masterpiece collection, the work of art known as '{work_of_art}'.")
    print("\nHere are the verses:\n")
    print(verses)

find_poet_and_work()