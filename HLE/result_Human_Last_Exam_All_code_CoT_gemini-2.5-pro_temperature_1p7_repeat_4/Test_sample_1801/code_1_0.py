def find_enron_laureates():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in the 2000 Enron "Ask Why" commercial.
    """
    # The two Nobel laureates in question
    laureate1 = "Paul Krugman"
    laureate2 = "Myron Scholes"
    
    # Details from the question
    company = "Enron"
    year = 2000
    tagline = "Ask Why"
    
    # Printing the final answer sentence
    print(f"The two Nobel laureates featured in an infamous {year} commercial for {company}, with the tagline \"{tagline}\", were:")
    print(f"- {laureate1}")
    print(f"- {laureate2}")

# Execute the function
find_enron_laureates()