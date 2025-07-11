def find_poet_and_work():
    """
    This function identifies and prints the author and the work of art
    associated with the provided poetic verses.
    """
    poet = "Rafael Alberti"
    poem_title = "'Ángel de los números' (Angel of the Numbers)"
    work_of_art = "the poetry collection 'Sobre los ángeles' (Concerning the Angels), published in 1929."

    print("The Spanish poet who wrote these verses is:")
    print(f"- {poet}")
    print("\nThe verses are from his poem titled:")
    print(f"- {poem_title}")
    print("\nThis poem is part of the larger work of art to which it was directed, which is:")
    print(f"- {work_of_art}")

if __name__ == "__main__":
    find_poet_and_work()