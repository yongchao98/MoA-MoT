def identify_coin_origin():
    """
    Identifies the origin of a French colonial coin based on its features.
    
    The coin is inscribed with "COLONIES FRANCOISES" and dates to 1721-1722.
    This specific issue of 9 Deniers was primarily minted for use in a
    particular North American colony to support the economic activities
    of the Company of the West.
    """
    
    coin_details = {
        "inscription": "COLONIES FRANCOISES",
        "era": "Louis XV (1721-1722)",
        "denomination": "9 Deniers",
        "primary_colony": "French Louisiana"
    }
    
    colony = coin_details["primary_colony"]
    
    print(f"This coin was primarily intended for the colony of: {colony}")

identify_coin_origin()