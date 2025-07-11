def find_enron_laureates():
    """
    Searches a predefined knowledge base to find the Nobel laureates
    featured in a specific Enron commercial.
    """
    # A small knowledge base representing information about commercials.
    commercial_database = {
        "Enron": [
            {
                "year": 2000,
                "tagline": "Ask Why",
                "featured_nobel_laureates": ["Paul Krugman", "Myron Scholes"]
            }
        ],
        "Apple": [
            {
                "year": 1997,
                "tagline": "Think Different",
                "featured_nobel_laureates": ["Albert Einstein", "Bob Dylan"] # Note: Einstein shown, Dylan a featured voice.
            }
        ]
    }

    company_to_find = "Enron"
    year_to_find = 2000
    tagline_to_find = "Ask Why"
    
    found_laureates = []
    
    # Search the database
    if company_to_find in commercial_database:
        for ad in commercial_database[company_to_find]:
            if ad.get("year") == year_to_find and ad.get("tagline") == tagline_to_find:
                found_laureates = ad.get("featured_nobel_laureates", [])
                break

    if found_laureates:
        print(f"The two Nobel laureates featured in the infamous {year_to_find} Enron commercial with the tagline '{tagline_to_find}' were:")
        for laureate in found_laureates:
            print(f"- {laureate}")
    else:
        print("Could not find the requested information in the database.")

find_enron_laureates()