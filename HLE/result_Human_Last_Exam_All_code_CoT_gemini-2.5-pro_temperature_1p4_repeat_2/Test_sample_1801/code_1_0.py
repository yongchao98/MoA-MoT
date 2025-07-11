def solve_enron_commercial_question():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in the 2000 Enron "Ask Why" commercial.
    """
    # Storing the names of the Nobel laureates
    laureate1 = "Paul Krugman"
    laureate2 = "Myron Scholes"

    # Storing the other details from the question for context
    company = "Enron"
    year = 2000
    tagline = "Ask Why"

    # Printing the final answer in a full sentence
    print(f"The two Nobel laureates featured in the infamous {year} commercial for {company} with the tagline '{tagline}' were:")
    print(f"- {laureate1}")
    print(f"- {laureate2}")

solve_enron_commercial_question()