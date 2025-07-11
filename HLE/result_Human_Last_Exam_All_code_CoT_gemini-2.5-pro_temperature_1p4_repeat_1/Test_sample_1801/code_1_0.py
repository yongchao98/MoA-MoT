def solve_enron_commercial_question():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in an infamous 2000 Enron commercial.
    """
    
    # The two individuals featured in the commercial.
    # Note: Paul Krugman was not a Nobel laureate at the time of the commercial in 2000,
    # but was awarded the prize in 2008. Merton Miller was a laureate, having won in 1990.
    laureate1 = "Paul Krugman"
    laureate2 = "Merton Miller"
    
    # Details from the user's question.
    company = "Enron"
    year = 2000
    tagline = "Ask Why"

    # Print the answer in a full sentence.
    print(f"The two Nobel laureates featured in the infamous {year} commercial for {company} with the tagline '{tagline}' were {laureate1} and {laureate2}.")

solve_enron_commercial_question()