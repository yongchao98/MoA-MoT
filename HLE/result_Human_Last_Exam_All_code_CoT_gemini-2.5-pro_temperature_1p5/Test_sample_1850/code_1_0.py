import re
import requests

def count_saints_in_paradise_lost():
    """
    Downloads Paradise Lost from Project Gutenberg and counts the number of
    uniquely named historical saints explicitly mentioned.
    """
    # The URL for the plain text version of Paradise Lost from Project Gutenberg.
    url = "https://www.gutenberg.org/cache/epub/20/pg20.txt"
    
    print(f"Attempting to download Paradise Lost from {url}...")
    
    try:
        # Fetch the content of the book
        response = requests.get(url, timeout=10)
        # Raise an exception if the download failed
        response.raise_for_status()
        text = response.text
        print("Download successful. Analyzing text...")

        # Use regular expression to find all instances of "Saint" followed by a capitalized name.
        # This is the most reliable way to find figures explicitly called "Saint".
        saints_found = re.findall(r"Saint\s+([A-Z][a-z]+)", text)
        
        # Get a unique set of names, as a saint might be mentioned more than once.
        unique_saints = sorted(list(set(saints_found)))
        
        count = len(unique_saints)
        
        print("\nFound the following historical saint(s) explicitly named:")
        if not unique_saints:
            print("None")
        else:
            for saint in unique_saints:
                print(f"- {saint}")

        print("\n---")
        # As requested, showing the final equation.
        print("Final Equation:")
        equation_str = " + ".join(["1"] * count) if count > 0 else "0"
        print(f"{equation_str} = {count}")
        print(f"Total number of historical saints mentioned by name: {count}")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. Please check your internet connection.")
        print(f"Details: {e}")

if __name__ == '__main__':
    count_saints_in_paradise_lost()