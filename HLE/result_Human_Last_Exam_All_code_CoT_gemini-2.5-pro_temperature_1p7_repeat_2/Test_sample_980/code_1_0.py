def calculate_station_value(station_name):
    """
    Calculates the value of a station name by summing the alphabetical
    positions of its letters (A=1, B=2, ...).
    """
    total_value = 0
    # Remove spaces and convert to uppercase to handle names like "Arnos Grove"
    cleaned_name = ''.join(filter(str.isalpha, station_name)).upper()
    
    # We will print the equation for clarity
    equation_parts = []
    
    for char in cleaned_name:
        value = ord(char) - ord('A') + 1
        total_value += value
        equation_parts.append(str(value))
    
    print(f"Calculating value for station: {station_name}")
    print(f"Letters: {' '.join(list(cleaned_name))}")
    print(f"Values:  {' + '.join(equation_parts)} = {total_value}")
    return total_value

# Last known station in the sequence to verify the rule
last_station = "Bounds Green"
last_value = calculate_station_value(last_station)
print("-" * 20)

# The next station in the sequence
next_station = "Arnos Grove"
next_value = calculate_station_value(next_station)

print("-" * 20)
print(f"The next number in the sequence is: {next_value}")