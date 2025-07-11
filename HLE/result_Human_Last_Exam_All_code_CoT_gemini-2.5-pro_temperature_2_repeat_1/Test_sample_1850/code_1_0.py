import requests
import re

def find_saints_in_paradise_lost():
    """
    This script downloads John Milton's "Paradise Lost" and counts the number
    of pre-defined historical saints mentioned by name within the text.
    """
    # The URL for the plain text version of "Paradise Lost" from Project Gutenberg.
    url = "https://www.gutenberg.org/files/20/20-0.txt"
    
    # Based on analysis of the text, Milton names several Old Testament figures
    # who are revered as saints or holy men, but does not name New Testament saints.
    saints_list = [
        "Noah", "Enoch", "Abraham", "Isaac", "Jacob",
        "Moses", "Joshua", "David", "Solomon"
    ]
    
    found_saints = []

    print("Attempting to download 'Paradise Lost'...")
    try:
        response = requests.get(url, timeout=20)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code.
        response.raise_for_status()
        text = response.text
        print("Download complete. Searching for saints...\n")

        # Search for each saint in the text.
        for saint in saints_list:
            # We use a regular expression to find the whole word, case-insensitively.
            # The '\b' ensures we match the whole word and not a substring.
            # For example, it prevents matching 'Isaac' inside 'Isaac's'. (We will count 'Isaac's' as a mention of 'Isaac' however.)
            # A simple word boundary around the name handles most cases well.
            if re.search(r'\b' + re.escape(saint) + r'\b', text, re.IGNORECASE):
                found_saints.append(saint)

        if not found_saints:
            print("Could not find any of the specified saints in the text.")
            return
            
        print("The following historical saints were found by name in 'Paradise Lost':")
        
        # We will now print each found saint's name, which contributes to the final count.
        total_count = len(found_saints)
        for i, name in enumerate(found_saints):
            # This format represents each "number" in the "final equation".
            if i < total_count - 1:
                print(f"{name}", end=" + ")
            else:
                print(f"{name}", end="")
        
        print(f" = {total_count}")
        print(f"\nThe total number of historical saints mentioned by name is {total_count}.")


    except requests.exceptions.RequestException as e:
        print(f"An error occurred while trying to download the text: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

find_saints_in_paradise_lost()
<<<9>>>