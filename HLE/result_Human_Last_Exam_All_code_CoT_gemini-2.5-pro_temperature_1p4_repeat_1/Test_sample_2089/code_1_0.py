def solve_film_trivia():
    """
    Identifies the film based on a specific, memorable scene.
    """
    film_title = "A Simple Plan"
    year = 1998
    character = "Jacob Mitchell (played by Billy Bob Thornton)"

    explanation = (
        f"The Oscar-nominated film is '{film_title}' ({year}).\n\n"
        f"The character is {character}. Overcome with guilt, he orchestrates his own death. "
        "When the local sheriff pulls over on a snowy road and offers him a ride, Jacob says 'Thank you', "
        "and then pulls an unloaded gun, forcing the sheriff to fatally shoot him. "
        "The vehicle was the sheriff's car, which is commonly misremembered as a bus."
    )
    
    print(explanation)

solve_film_trivia()