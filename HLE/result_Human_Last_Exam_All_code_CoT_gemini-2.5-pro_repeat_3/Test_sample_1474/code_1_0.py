def solve_movie_mystery():
    """
    Analyzes the works of directors Fritz Lang and William Friedkin to find a common element.
    """
    print("Investigating common themes in the films of Fritz Lang and William Friedkin...")
    print("-" * 70)

    # Evaluation of Choice A: Aboriginal masks
    print("Choice A: Aboriginal masks")
    print("Analysis: This is not a known element in either director's major works. Friedkin's 'The Exorcist' features an Assyrian amulet, which is incorrect for this choice. Lang's films do not feature them.")
    print("Result: Incorrect.\n")

    # Evaluation of Choice B: Magic wands
    print("Choice B: Magic wands")
    print("Analysis: This classic fantasy trope does not fit the style of either director, who are known for noir, thrillers, and horror, not high fantasy.")
    print("Result: Incorrect.\n")

    # Evaluation of Choice C: The first ever cyborgs on screen
    print("Choice C: The first ever cyborgs on screen")
    print("Analysis: Fritz Lang's 'Metropolis' (1927) features the Maschinenmensch, a famous early android. However, William Friedkin's filmography does not contain cyborgs.")
    print("Result: Incorrect, as it does not apply to both directors.\n")

    # Evaluation of Choice D: Bugs
    print("Choice D: Bugs")
    print("Analysis: This theme appears in the work of both directors, both literally and metaphorically.")
    print("  - William Friedkin: Directed 'Bug' (2006), a film about paranoia and delusion involving insects. He also used flies as a motif for demonic presence in 'The Exorcist' (1973).")
    print("  - Fritz Lang: His final film, 'The Thousand Eyes of Dr. Mabuse' (1960), revolves around a network of spies using hidden cameras and microphones—technological 'bugs'.")
    print("The theme of 'bugs' connects both directors through imagery of paranoia, infestation, and hidden surveillance.")
    print("Result: Correct.\n")

    print("-" * 70)
    print("Final Conclusion Equation:")
    lang_movie_year = 1960
    friedkin_movie_year_1 = 1973
    friedkin_movie_year_2 = 2006
    print(f"Lang ('The Thousand Eyes of Dr. Mabuse', {lang_movie_year}) + Friedkin ('The Exorcist', {friedkin_movie_year_1}; 'Bug', {friedkin_movie_year_2})")
    print("=> Common Theme: Bugs (Metaphorical and Literal)")
    print("-" * 70)

solve_movie_mystery()