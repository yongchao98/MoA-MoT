def calculate_word_value(word):
    """Calculates the value of a word based on A=1, B=2, etc."""
    word = word.lower()
    value = 0
    # Store letter values for the final printout
    letter_values = []
    for char in word:
        if 'a' <= char <= 'z':
            char_val = ord(char) - ord('a') + 1
            value += char_val
            letter_values.append(str(char_val))
    return value, letter_values

def find_next_in_sequence():
    """Finds the 34th number in the sequence of US state capital values."""
    # List of states and capitals, sorted alphabetically by state
    states_capitals = [
        ("Alabama", "Montgomery"), ("Alaska", "Juneau"), ("Arizona", "Phoenix"),
        ("Arkansas", "Little Rock"), ("California", "Sacramento"), ("Colorado", "Denver"),
        ("Connecticut", "Hartford"), ("Delaware", "Dover"), ("Florida", "Tallahassee"),
        ("Georgia", "Atlanta"), ("Hawaii", "Honolulu"), ("Idaho", "Boise"),
        ("Illinois", "Springfield"), ("Indiana", "Indianapolis"), ("Iowa", "Des Moines"),
        ("Kansas", "Topeka"), ("Kentucky", "Frankfort"), ("Louisiana", "Baton Rouge"),
        ("Maine", "Augusta"), ("Maryland", "Annapolis"), ("Massachusetts", "Boston"),
        ("Michigan", "Lansing"), ("Minnesota", "Saint Paul"), ("Mississippi", "Jackson"),
        ("Missouri", "Jefferson City"), ("Montana", "Helena"), ("Nebraska", "Lincoln"),
        ("Nevada", "Carson City"), ("New Hampshire", "Concord"), ("New Jersey", "Trenton"),
        ("New Mexico", "Santa Fe"), ("New York", "Albany"), ("North Carolina", "Raleigh"),
        ("North Dakota", "Bismarck"), ("Ohio", "Columbus"), ("Oklahoma", "Oklahoma City"),
        ("Oregon", "Salem"), ("Pennsylvania", "Harrisburg"), ("Rhode Island", "Providence"),
        ("South Carolina", "Columbia"), ("South Dakota", "Pierre"), ("Tennessee", "Nashville"),
        ("Texas", "Austin"), ("Utah", "Salt Lake City"), ("Vermont", "Montpelier"),
        ("Virginia", "Richmond"), ("Washington", "Olympia"), ("West Virginia", "Charleston"),
        ("Wisconsin", "Madison"), ("Wyoming", "Cheyenne")
    ]

    # The sequence has 33 numbers, so we need the 34th.
    # The 34th state in the alphabetical list is Ohio, capital Columbus.
    target_state, target_capital = states_capitals[33]
    
    # Handle two-word capitals correctly by removing spaces
    capital_to_process = target_capital.replace(" ", "")

    value, letter_values = calculate_word_value(capital_to_process)
    
    # Format the output equation
    letters = " + ".join(list(capital_to_process.upper()))
    values_str = " + ".join(letter_values)
    
    print(f"The next number is derived from the capital of the 34th state, {target_state}: {target_capital}.")
    print(f"Calculation for '{capital_to_process.upper()}':")
    print(f"{letters}")
    print(f"{values_str} = {value}")

find_next_in_sequence()