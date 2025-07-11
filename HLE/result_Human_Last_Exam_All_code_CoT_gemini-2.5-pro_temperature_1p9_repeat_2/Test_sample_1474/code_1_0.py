def find_common_theme():
    """
    Analyzes the oeuvres of Fritz Lang and William Friedkin to identify a common theme.
    """

    print("Investigating a common theme in the films of Fritz Lang and William Friedkin...")
    print("="*60)

    # Step 1: Analyze Fritz Lang's work for the theme of 'Bugs'
    lang_film = "Die Spinnen (The Spiders)"
    lang_year = 1919
    print(f"Director: Fritz Lang")
    print(f"Evidence: Lang directed '{lang_film}' in {lang_year}.")
    print("This film's title translates to 'The Spiders'. Spiders, while technically arachnids, are colloquially known as 'bugs'.")
    print("-"*60)

    # Step 2: Analyze William Friedkin's work for the theme of 'Bugs'
    friedkin_film_1 = "The Exorcist"
    friedkin_year_1 = 1973
    friedkin_film_2 = "Bug"
    friedkin_year_2 = 2006
    print(f"Director: William Friedkin")
    print(f"Evidence 1: In '{friedkin_film_1}' ({friedkin_year_1}), the demon Pazuzu is associated with bringing locust swarms, which are bugs.")
    print(f"Evidence 2: Friedkin also directed the film '{friedkin_film_2}' ({friedkin_year_2}), which centers on a paranoid delusion of a bug infestation.")
    print("="*60)
    
    # Step 3: Conclude and present the 'equation' of evidence based on film years
    print("Conclusion: The theme of 'Bugs' is present in the work of both directors.")
    print("\nThe numerical evidence, based on the years of the relevant films, is as follows:")
    
    # This fulfills the prompt's requirement to output each number in the final equation.
    print(f"Fritz Lang's evidence year: {lang_year}")
    print(f"William Friedkin's evidence year 1: {friedkin_year_1}")
    print(f"William Friedkin's evidence year 2: {friedkin_year_2}")

find_common_theme()