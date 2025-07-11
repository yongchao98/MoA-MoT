def solve_musical_joke():
    """
    This function provides the details of a famous musical joke in a piano concerto.
    The information is based on established musicological analysis.
    """
    composer_surname = "Beethoven"
    opus_number = 73
    
    # The flourish begins with a 'false entry' in measure 97.
    start_measure = 97
    
    # The cellos recover by skipping a beat (resting) in measure 100.
    end_measure = 100
    
    # Format the measure range as a hyphenated string
    measure_range = f"{start_measure}-{end_measure}"
    
    # Combine the elements into the final comma-separated string
    # The problem asks to output each number in the final equation. 
    # Here, "equation" is interpreted as the final formatted string.
    print(f"{composer_surname}, {opus_number}, {measure_range}")

solve_musical_joke()