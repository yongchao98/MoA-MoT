def get_historical_judge():
    """
    This function identifies and prints the name of the judge of Bártfa (now Bardejov) in 1461.
    The information is based on historical data.
    """
    # Define the known historical data
    year = 1461
    city_name_hungarian = "Bártfa"
    city_name_slovak = "Bardejov"
    judge_name = "Alexius"

    # The prompt requires outputting the number, so we will include the year in the final sentence.
    # The equation here is the statement of fact.
    print(f"In the year {year}, the judge of {city_name_hungarian} (now known as {city_name_slovak}) was {judge_name}.")

# Call the function to display the answer
get_historical_judge()