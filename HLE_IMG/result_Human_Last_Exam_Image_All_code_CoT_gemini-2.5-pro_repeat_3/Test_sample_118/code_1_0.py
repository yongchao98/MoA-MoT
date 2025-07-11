def find_best_path():
    """
    This function determines the best path from Billoweg to Laubiweg
    based on the provided criteria and prints the result in the specified format.
    The analysis and station counts are pre-determined by examining the map.
    """

    # Path components based on the step-by-step analysis
    start_symbol = "S"
    end_symbol = "E"
    
    first_line_number = 72
    exchange_station_name = "Guggachstr"
    second_line_number = 9
    total_station_count = 23

    # Construct the final output string as per the specified format.
    # The numbers are explicitly part of the output string's construction.
    result_path = f"{start_symbol} - {first_line_number} - {exchange_station_name} - {second_line_number} - {end_symbol}; {total_station_count}"

    print(result_path)

find_best_path()