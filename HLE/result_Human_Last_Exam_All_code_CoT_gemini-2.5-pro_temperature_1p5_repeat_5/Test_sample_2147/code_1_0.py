def calculate_age_at_event():
    """
    This function calculates and prints Pavel Milyukov's age
    at the start of the "World of Art" movement.
    """
    milyukov_birth_year = 1859
    world_of_art_start_year = 1898

    # Calculate Milyukov's age at the time of the first exhibition
    age_at_event = world_of_art_start_year - milyukov_birth_year

    print("Calculating Pavel Milyukov's age when the 'World of Art' movement began in 1898.")
    # The prompt requires printing each number in the final equation.
    print(f"Equation: {world_of_art_start_year} - {milyukov_birth_year} = {age_at_event}")
    print(f"Pavel Milyukov was {age_at_event} years old.")

calculate_age_at_event()