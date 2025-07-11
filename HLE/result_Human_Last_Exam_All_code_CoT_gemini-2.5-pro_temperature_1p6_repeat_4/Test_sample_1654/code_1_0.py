import requests
import re

def count_coat_changes():
    """
    Downloads and analyzes "The Chameleon" by Anton Chekov to count
    the symbolic descriptions of Otchumyelov's coat changing.
    """
    try:
        # Step 1: Download the story from a public domain source (Project Gutenberg).
        # This plain text file contains the standard Constance Garnett translation.
        url = "https://www.gutenberg.org/files/13415/13415-0.txt"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text

        # Step 2: Define regular expression patterns for taking the coat ON and OFF.
        # These phrases are key indicators of Otchumyelov's changing state of mind.
        
        # Pattern for taking the coat OFF (he feels "hot"). The phrase "off with my coat" appears twice.
        off_pattern = r"off with my coat"
        
        # Pattern for putting the coat ON (he feels "chilly" or reasserts authority).
        # We combine two distinct phrases into one pattern using '|' (OR).
        on_pattern = r"(Put on my coat|wraps himself in his greatcoat)"

        # Step 3: Find all occurrences of these patterns in the text.
        off_matches = re.findall(off_pattern, text, re.IGNORECASE)
        on_matches = re.findall(on_pattern, text, re.IGNORECASE)

        num_off = len(off_matches)
        num_on = len(on_matches)
        total_changes = num_off + num_on

        # Step 4: Print the results and the final equation.
        print(f"Analysis of Otchumyelov's coat in 'The Chameleon':\n")
        print(f"Found {num_off} instance(s) of the coat being taken OFF.")
        print(f"Found {num_on} instance(s) of the coat being put ON.")
        print("-" * 40)
        print("The total number of symbolic descriptions is calculated as follows:")
        print(f"{num_off} + {num_on} = {total_changes}")
        print("-" * 40)

    except requests.exceptions.RequestException as e:
        print(f"Failed to download the story. Error: {e}")
        print("Cannot perform analysis. Please check your internet connection.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    count_coat_changes()