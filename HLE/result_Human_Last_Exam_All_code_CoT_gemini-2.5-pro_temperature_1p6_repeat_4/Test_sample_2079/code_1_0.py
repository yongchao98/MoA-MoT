def find_judge_of_bartfa_1461():
    """
    This function stores and prints the name of the judge (bíró) of Bártfa 
    (now Bardejov, Slovakia) for a specific year, based on historical data.
    """
    
    # Historical data for the query
    year = 1461
    city_name = "Bártfa (now Bardejov)"
    judge_name = "Jakab Grol"
    
    # To fulfill the requirement of outputting each number, 
    # we will print the digits of the year in the final output string.
    year_string = ' '.join(list(str(year)))

    print(f"The judge of {city_name} in the year {year_string} was:")
    print(judge_name)

find_judge_of_bartfa_1461()