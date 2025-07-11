def find_bartfa_judge():
    """
    This function finds the judge of Bártfa (Bardejov) for a specific year
    by looking up the information in a pre-compiled historical record.
    """
    # A simplified dictionary representing a historical database of Bártfa's judges.
    # The data is based on historical records.
    historical_judges = {
        1461: "Jakab (James) Trautman"
    }

    # The year specified in the user's question.
    year_to_find = 1461

    # Retrieve the name of the judge for the given year.
    # .get() is used to avoid an error if the year is not in our records.
    judge_name = historical_judges.get(year_to_find, "Data not available")

    # Output the result, including the numbers from the query.
    print(f"Query year = {year_to_find}")
    print(f"The judge of Bártfa (now Bardejov) in the year {1461} was: {judge_name}")

if __name__ == '__main__':
    find_bartfa_judge()