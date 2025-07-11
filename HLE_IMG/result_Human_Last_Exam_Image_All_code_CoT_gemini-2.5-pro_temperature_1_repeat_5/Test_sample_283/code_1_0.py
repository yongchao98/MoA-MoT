def solve_historical_trends():
    """
    This function analyzes the provided Google Trends data to identify the country
    where three specific historical events are commemorated.
    """
    
    # Analysis based on the timing of recurring peaks in the graph
    event_red_line = "Martin Luther King Jr. Day (celebrated in late January)"
    event_yellow_line = "Juneteenth (celebrated on June 19th)"
    event_blue_line = "Memorial Day (celebrated on the last Monday of May)"
    
    country = "United States"
    
    print("Analysis of the Google Trends Graph:")
    print(f"1. The red line shows annual peaks in late January, which corresponds to '{event_red_line}'.")
    print(f"2. The yellow line shows sharp annual peaks in mid-June, corresponding to '{event_yellow_line}'.")
    print(f"3. The blue line shows annual peaks at the end of May, corresponding to '{event_blue_line}'.")
    
    print("\nConclusion:")
    print(f"These three major historical commemorations are all federal holidays in the {country}.")

solve_historical_trends()
