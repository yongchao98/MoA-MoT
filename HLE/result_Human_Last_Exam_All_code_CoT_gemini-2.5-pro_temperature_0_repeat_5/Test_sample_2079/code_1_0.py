def find_judge_of_bartfa():
    """
    This function provides the name of the judge of Bártfa (Bardejov) in 1461.
    """
    # Information based on historical records for the specified year.
    city_name_historical = "Bártfa"
    city_name_modern = "Bardejov"
    year = 1461
    judge_name = "Leonard of Uničov"

    # The final output string, constructed from the variables.
    # The number from the query is explicitly included in the output.
    answer = f"In the year {year}, the judge of {city_name_historical} (now {city_name_modern}) was {judge_name}."

    print(answer)

# Execute the function to print the answer.
find_judge_of_bartfa()