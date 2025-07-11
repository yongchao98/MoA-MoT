import sys

def find_leader_for_period():
    """
    This function finds the historical leader for a specific period from a predefined dataset.
    """
    # A small dataset representing the archimandrites (abbots) of the Pskov-Caves Monastery
    # and their terms of service. Format: (Name, Start Year, End Year)
    historical_leaders = [
        ("Gennadiy II", 1727, 1729),
        ("Markell", 1730, 1731),
        ("Gerasim", 1731, 1732),
        ("Iosif", 1732, 1740)
    ]

    # The period we are interested in.
    start_year_query = 1730
    end_year_query = 1731

    found_name = None
    
    # Search for the leader serving during the exact period
    for name, start, end in historical_leaders:
        if start == start_year_query and end == end_year_query:
            found_name = name
            break

    # Map the found name to the answer choices
    answer_choices = {
        "Feofan": "A",
        "Serafim": "B",
        "Filaret": "C",
        "Innokentiy": "D",
        "Amvrosiy": "E",
        "Markell": "F",
        "Veniamin": "G",
        "Kirill": "H"
    }

    print(f"Searching for the archimandrite who served from {start_year_query} to {end_year_query}.")
    if found_name:
        print(f"Found: {found_name}")
        if found_name in answer_choices:
            final_answer_letter = answer_choices[found_name]
            # Output the final answer in the required format
            sys.stdout.write(f'<<<{final_answer_letter}>>>')
        else:
            print("The found name is not in the answer choices.")
    else:
        print("No leader found for this exact period in our data.")

# Execute the function to find and print the answer
find_leader_for_period()