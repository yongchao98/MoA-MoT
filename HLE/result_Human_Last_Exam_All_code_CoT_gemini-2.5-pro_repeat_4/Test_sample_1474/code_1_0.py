def find_common_element():
    """
    Analyzes the works of Fritz Lang and William Friedkin to find a shared element.
    """
    print("Finding a shared cinematic element in the works of Fritz Lang and William Friedkin...")
    print("-" * 75)

    # Step 1: Examine Fritz Lang's filmography for the element 'Bugs'.
    lang_film = "Dr. Mabuse the Gambler"
    lang_year = 1922
    lang_context = "a character under hypnosis hallucinates a terrifying swarm of insects."
    print(f"1. In Fritz Lang's '{lang_film}', from the year {lang_year}, {lang_context}")

    # Step 2: Examine William Friedkin's filmography for the element 'Bugs'.
    friedkin_film = "Bug"
    friedkin_year = 2006
    friedkin_context = "the entire plot revolves around a paranoia of a government-implanted insect infestation."
    print(f"2. In William Friedkin's '{friedkin_film}', from the year {friedkin_year}, {friedkin_context}")
    print("-" * 75)

    # Step 3: Conclude that 'Bugs' is the common element.
    print("Conclusion: Both directors used 'bugs' as a powerful motif for psychological horror and paranoia.")
    print("\nThis can be represented as the following 'equation':")

    # Step 4: Display the final answer in the requested format.
    print(f"Use of bugs by Lang in {lang_year} + Use of bugs by Friedkin in {friedkin_year} = Common Theme: 'Bugs'")

find_common_element()