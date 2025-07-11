def find_enron_laureates():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in Enron's 2000 "Ask Why" commercial.
    """
    # The two Nobel laureates featured in the commercial
    laureates = ["Paul Krugman", "Myron Scholes"]
    
    # Details of the commercial
    company = "Enron"
    tagline = "Ask Why"
    year = 2000
    
    # Print the final answer in a full sentence
    print(f"The two Nobel laureates featured in the infamous {year} commercial for {company}, with the tagline '{tagline}', were {laureates[0]} and {laureates[1]}.")

find_enron_laureates()