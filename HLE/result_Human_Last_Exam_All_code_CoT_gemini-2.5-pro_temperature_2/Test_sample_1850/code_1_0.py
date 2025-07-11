import requests
import re

def count_saints_in_paradise_lost():
    """
    Downloads Milton's Paradise Lost and counts the number of historical saints
    mentioned by name.
    """
    # URL for the plain text of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/files/20/20-0.txt"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text: {e}")
        return

    # A list of historical figures considered saints who are mentioned by name in the poem.
    # This list is based on literary analysis of the text.
    # It includes Old Testament figures, New Testament figures, and post-biblical figures.
    saints_to_find = [
        "Adam", "Eve", "Enoch", "Abraham", "Moses", "Joshua",
        "David", "Solomon", "John", "Mary", "Dominic", "Francis"
    ]
    
    found_saints = []
    
    # Search for each saint's name in the text
    for saint in saints_to_find:
        # We use word boundaries (\b) to ensure we match the whole name
        # e.g., to not match "Adam" in "madam".
        if re.search(r'\b' + re.escape(saint) + r'\b', text, re.IGNORECASE):
            found_saints.append(saint)

    if not found_saints:
        print("No historical saints from the list were found in the text.")
        return

    # Sort the list for consistent output
    found_saints.sort()
    
    # Print the equation as requested
    equation = " + ".join(found_saints)
    total = len(found_saints)
    
    print("Based on analysis of the text, the following historical saints are mentioned by name:")
    print(f"\n{equation} = {total}\n")
    print(f"The total number of historical saints mentioned by name is {total}.")


if __name__ == '__main__':
    count_saints_in_paradise_lost()