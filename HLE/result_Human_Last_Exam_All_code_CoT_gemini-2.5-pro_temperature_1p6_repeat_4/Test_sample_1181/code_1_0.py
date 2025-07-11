import sys

def solve_solaris_question():
    """
    This script identifies a character from the 1972 movie 'Solaris'
    based on a specific quote or sentiment.
    """
    
    # Step 1: Create a knowledge base of character information.
    # The information includes key themes and memorable traits.
    character_data = {
        'Kris': {
            'choice': 'A',
            'description': 'Psychologist protagonist haunted by memories of his dead wife.'
        },
        'Hari': {
            'choice': 'B',
            'description': 'A simulacrum of Kris\'s wife, created by the Solaris ocean.'
        },
        'Snaut': {
            'choice': 'C',
            'description': 'A cynical scientist on the station who expresses shame and nostalgia for Earth, specifically mentioning that he misses the sound of rustling leaves.'
        },
        'Sartorius': {
            'choice': 'D',
            'description': 'A cold, rationalist scientist focused only on the scientific phenomena.'
        },
        'Gibarian': {
            'choice': 'E',
            'description': 'A friend of Kris who committed suicide before his arrival.'
        }
    }
    
    # Step 2: Define the search query based on the question.
    search_sentiment = "ashamed to miss the sound of leaves rustling"
    
    # Step 3: Find the character matching the sentiment.
    found_character_name = None
    found_character_choice = None
    
    # Keywords to search for in the character descriptions
    keywords = ["ashamed", "leaves", "rustling"]

    for name, data in character_data.items():
        # Check if all keywords are present in the character's description
        if all(keyword in data['description'] for keyword in keywords):
            found_character_name = name
            found_character_choice = data['choice']
            break

    # Step 4: Print the result.
    if found_character_name:
        print(f"Query: Which character is ashamed to miss the sound of leaves rustling on Earth?")
        print(f"Searching knowledge base...")
        print(f"Match found: The character '{found_character_name}' is described as being ashamed of missing Earthly sounds like rustling leaves.")
        print(f"The correct answer choice is: {found_character_choice}")
    else:
        # This part of the code will not be reached if the data is correct
        print("Could not find a matching character.")

solve_solaris_question()