import requests
import re

def count_coat_mentions():
    """
    Fetches Anton Chekov's "The Chameleon" and counts the mentions
    of Otchumyelov's coat.
    """
    try:
        # URL for "The Wife, and Other Stories" by Anton Chekhov (trans. Constance Garnett)
        # from Project Gutenberg, which includes "The Chameleon"
        url = "https://www.gutenberg.org/cache/epub/13509/pg13509.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        full_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching the story: {e}")
        return

    try:
        # Step 1: Isolate the text of "The Chameleon" from the collection.
        # The story starts with its title and ends just before the next story's title.
        start_marker = "THE CHAMELEON"
        # The next story in this specific collection is "BIG VOLODYA AND LITTLE VOLODYA"
        end_marker = "BIG VOLODYA AND LITTLE VOLODYA"
        
        start_index = full_text.index(start_marker)
        end_index = full_text.index(end_marker, start_index)
        
        story_text = full_text[start_index:end_index]
    except ValueError:
        print("Error: Could not find the start or end of the story in the downloaded text.")
        return

    # Step 2: Convert text to lowercase for case-insensitive counting.
    story_text_lower = story_text.lower()
    
    # Step 3: Count occurrences of "great-coat" and standalone "coat".
    # Using word boundaries (\b) ensures we match whole words and that counting
    # "\bcoat\b" does not also count "great-coat".
    great_coat_count = len(re.findall(r'\bgreat-coat\b', story_text_lower))
    coat_count = len(re.findall(r'\bcoat\b', story_text_lower))
    
    total_count = great_coat_count + coat_count
    
    # Step 4: Output the results, showing each number in the final equation.
    print("Analysis of the word 'coat' in Chekov's 'The Chameleon' (Constance Garnett translation):")
    print(f"Number of times 'great-coat' is mentioned: {great_coat_count}")
    print(f"Number of times 'coat' is mentioned (standalone): {coat_count}")
    print(f"Total descriptions of the coat: {great_coat_count} + {coat_count} = {total_count}")

if __name__ == '__main__':
    count_coat_mentions()