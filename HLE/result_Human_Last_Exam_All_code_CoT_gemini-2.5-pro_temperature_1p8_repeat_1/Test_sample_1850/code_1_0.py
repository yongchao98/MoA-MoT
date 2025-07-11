import requests
import re

def count_historical_saints_in_paradise_lost():
    """
    This script searches for the names of historical, post-biblical saints
    in the text of Milton's 'Paradise Lost' to determine how many are mentioned.
    """
    print("This script will check for historical (post-biblical) saints in 'Paradise Lost'.")
    print("Step 1: Downloading the text of 'Paradise Lost'...")

    # The URL for the plain text of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/cache/epub/20/pg20.txt"
    text = ""
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()  # This will raise an HTTPError for bad responses
        text = response.text
        print("Successfully downloaded 'Paradise Lost'.")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve the text from the internet. {e}")
        print("Cannot proceed without the text.")
        return

    # A representative list of major historical saints who lived after the biblical era.
    # Note: Biblical figures like Peter, David, Moses etc., are excluded as the focus
    # is on saints from post-apostolic church history.
    saints_to_check = [
        "Augustine", "Ambrose", "Jerome", "Gregory", "Benedict", "Francis",
        "Dominic", "Aquinas", "Catherine", "Teresa", "Ignatius", "George", "Patrick"
    ]

    print("\nStep 2: Searching the text for the following names:")
    print(", ".join(saints_to_check))
    print("-" * 50)

    found_saints = []
    text_lower = text.lower()

    for saint in saints_to_check:
        # Use regex with word boundaries (\b) to match only whole words, case-insensitively.
        # This prevents matching "francis" in "francs" or "dominic" in "dominion".
        pattern = r'\b' + re.escape(saint.lower()) + r'\b'
        if re.search(pattern, text_lower):
            found_saints.append(saint)

    # Step 3: Report the findings and present the final equation.
    print("Step 3: Compiling results...")
    total_found = len(found_saints)

    if total_found == 0:
        print("\nConclusion: As expected from Milton's theology and the poem's subject,")
        print("no historical, post-biblical saints from the list were mentioned by name.")
    else:
        print(f"\nFound {total_found} historical saint(s) mentioned by name: {', '.join(found_saints)}")

    # Fulfilling the request to output each number in the final equation.
    # If saints were found, this would be a sum of 1s.
    # Since none are found, the set of numbers is empty, and the total is 0.
    print("\nFinal Equation:")
    if total_found > 0:
        # This part of the code is unlikely to run, but is included for completeness.
        equation_numbers = ["1"] * total_found
        equation_string = " + ".join(equation_numbers)
        print(f"{equation_string} = {total_found}")
    else:
        # If no saints are found, the equation is simply a statement of the total.
        print(f"0 = {total_found}")


if __name__ == '__main__':
    count_historical_saints_in_paradise_lost()

<<<0>>>