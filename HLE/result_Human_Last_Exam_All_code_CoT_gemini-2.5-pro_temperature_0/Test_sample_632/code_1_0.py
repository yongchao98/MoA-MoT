def find_historical_feature():
    """
    This function provides information about a feature removed from
    Piazza della Rotonda around 1950.
    """
    year_of_removal = 1950
    feature_name = "The tram lines"
    location = "Piazza della Rotonda, Rome"
    
    print(f"The unique feature removed from the {location} around the year {year_of_removal} was:")
    print(f">>> {feature_name}")
    
    explanation = (
        "For the first half of the 20th century, tram lines ran directly through the piazza, "
        "passing in front of the Pantheon. These lines were a significant part of the city's landscape. "
        "As part of a post-war modernization plan and to reduce vibrations and pollution near the ancient monument, "
        "the city administration removed the tramway from the historic center, including the tracks in the piazza."
    )
    
    print("\nExplanation:")
    print(explanation)

# Execute the function to display the answer.
find_historical_feature()