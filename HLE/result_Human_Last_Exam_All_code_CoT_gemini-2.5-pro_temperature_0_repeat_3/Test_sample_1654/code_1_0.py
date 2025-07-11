import requests
import re

def count_coat_descriptions():
    """
    This script finds the number of times Chekov described Otchumyelov's coat
    in the short story "The Chameleon".
    """
    try:
        # Step 1: Fetch the collection of short stories from Project Gutenberg.
        # This file contains "The Chameleon".
        url = "https://www.gutenberg.org/files/13544/13544-0.txt"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes (like 404, 500).
        full_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not fetch the story text from the internet. Please check your connection.")
        print(f"Details: {e}")
        return

    try:
        # Step 2: Isolate the text of "The Chameleon".
        # The story begins after its title heading and ends before the next story's title.
        start_marker = "THE CHAMELEON"
        end_marker = "THE BEGGAR"
        
        # Find the start of the story using a regular expression to match the title as a whole line.
        start_match = re.search(r'^\s*' + re.escape(start_marker) + r'\s*$', full_text, re.MULTILINE)
        if not start_match:
            print(f"Error: Could not find the start of the story '{start_marker}'.")
            return
        
        start_index = start_match.end()
        
        # Find the end of the story, which is the beginning of the next one.
        end_index = full_text.find(end_marker, start_index)
        if end_index == -1:
            print(f"Error: Could not find the end of the story '{start_marker}'.")
            return
            
        story_text = full_text[start_index:end_index]
        
    except Exception as e:
        print(f"An unexpected error occurred while parsing the text: {e}")
        return

    # Step 3: Count the specific descriptions of the coat or symbolic actions involving it.
    # We convert the text to lowercase to ensure the search is case-insensitive.
    lower_story_text = story_text.lower()

    # Phrase 1: The initial description of the coat.
    desc1_phrase = "wearing a new overcoat"
    count1 = lower_story_text.count(desc1_phrase)

    # Phrase 2: Taking the coat off when he thinks the dog belongs to the General.
    desc2_phrase = "help me off with my coat"
    count2 = lower_story_text.count(desc2_phrase)

    # Phrase 3: Putting the coat back on when he thinks the dog is a stray.
    desc3_phrase = "help me on with my coat"
    count3 = lower_story_text.count(desc3_phrase)

    # Phrase 4: The final description of him wrapping himself in his "greatcoat".
    desc4_phrase = "wrapping himself in his greatcoat"
    count4 = lower_story_text.count(desc4_phrase)

    total_count = count1 + count2 + count3 + count4

    # Step 4: Print the results, including the final equation as requested.
    print("Counting the symbolic descriptions of Otchumyelov's coat in 'The Chameleon':")
    print(f"\n1. The initial description ('...{desc1_phrase}...'): {count1} time")
    print(f"2. Taking the coat off ('...{desc2_phrase}...'): {count2} time")
    print(f"3. Putting the coat on ('...{desc3_phrase}...'): {count3} time")
    print(f"4. The final wrapping action ('...{desc4_phrase}...'): {count4} time")
    
    print("\nThe final equation is:")
    print(f"{count1} + {count2} + {count3} + {count4} = {total_count}")
    
    print(f"\nTherefore, Chekov described Otchumyelov's coat {total_count} times.")

# Execute the function
count_coat_descriptions()