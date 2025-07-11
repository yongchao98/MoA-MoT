def solve_enron_commercial_question():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in a 2000 Enron commercial.
    """
    # Details of the commercial
    company = "Enron"
    year = 2000
    tagline = "Ask Why"

    # The two Nobel laureates
    laureate_1 = "Paul Krugman"
    laureate_2 = "Merton Miller"

    # Print the answer in a structured way
    print(f"The two Nobel laureates featured in the infamous {year} commercial for {company} were:")
    print(f"1. {laureate_1}")
    print(f"2. {laureate_2}")
    print(f"The tagline for the commercial was \"{tagline}\".")

solve_enron_commercial_question()