def solve_trivia():
    """
    This function provides the answer to the historical trivia question
    about the Piazza della Rotonda.
    """
    # The unique feature that was removed by the city administration.
    removed_feature = "the tram lines"
    
    # The approximate year of removal.
    year_of_removal = 1950
    
    # The location in question.
    location = "Piazza della Rotonda"

    # Construct the answer string.
    answer = (f"The unique feature in the '{location}' that was removed around the year {year_of_removal} "
              f"was {removed_feature}.")

    print(answer)

solve_trivia()