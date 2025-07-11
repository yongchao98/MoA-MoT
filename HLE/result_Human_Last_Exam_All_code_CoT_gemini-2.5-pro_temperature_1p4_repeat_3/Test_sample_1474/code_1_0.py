def find_common_element():
    """
    Analyzes the filmographies of Fritz Lang and William Friedkin to find a shared thematic element.
    """
    director_1 = "Fritz Lang"
    director_2 = "William Friedkin"
    common_theme = "Bugs"

    # Evidence for Director 1
    lang_film = "Die Spinnen (The Spiders)"
    lang_connection = "Features a villainous criminal organization named 'The Spiders'."

    # Evidence for Director 2
    friedkin_film_1 = "The Exorcist"
    friedkin_film_2 = "Bug"
    friedkin_connection = f"In '{friedkin_film_1}', the demon Pazuzu is associated with locusts. '{friedkin_film_2}' is a thriller about insect-based paranoia."

    print("Step 1: Analyzing the works of Fritz Lang.")
    print(f"   - In the film '{lang_film}', the antagonists are a group named after an arachnid.")

    print("\nStep 2: Analyzing the works of William Friedkin.")
    print(f"   - In the film '{friedkin_film_1}', the demonic entity is linked to swarms of locusts.")
    print(f"   - The film '{friedkin_film_2}''s plot is driven by a delusion about a bug infestation.")

    print("\nStep 3: Forming the final conclusion.")
    print("   - Both directors use bugs/arachnids as a source of paranoia, evil, or chaos.")
    print("\nFinal Equation:")
    print(f"Lang's '{lang_film}' + Friedkin's '{friedkin_film_1}' ==> Common theme of '{common_theme}'")

find_common_element()