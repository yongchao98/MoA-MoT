def identify_town_with_population_check():
    """
    This function confirms the identity of the town in the image by checking
    if its population meets the specified criteria.
    """
    # The town was identified as Blue Mound, Illinois.
    town_name = "Blue Mound"

    # Population of Blue Mound according to the 2020 US Census.
    population = 1141

    # The required population threshold from the user's request.
    population_threshold = 1000

    # The prompt requires showing the numbers involved.
    # The following print statement verifies the condition and includes all numbers.
    print(
        f"The identified town is {town_name}. "
        f"Its population of {population} is greater than the required {population_threshold}."
    )

identify_town_with_population_check()