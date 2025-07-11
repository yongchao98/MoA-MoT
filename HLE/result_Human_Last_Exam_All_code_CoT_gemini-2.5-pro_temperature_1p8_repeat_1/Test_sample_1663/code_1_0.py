def solve_poetry_query():
    """
    Identifies the author and work for the provided Spanish poetry verses.
    """
    poet = "Rafael Alberti"
    work_of_art = "Sobre los ángeles (Concerning the Angels)"
    poem_title = "El ángel de los números (The Angel of Numbers)"
    year_published = 1929

    print("The Spanish poet who wrote these verses is:")
    print(f"- {poet}")
    print("\nThe verses are from the poem '" + poem_title + "'.")
    print(f"\nThis poem is part of the work of art (a collection of poems) titled:")
    print(f"- {work_of_art}")
    print("\nThis seminal work of Spanish Surrealism was published in the year:")
    print(f"- {year_published}")

solve_poetry_query()