import requests
import re

def count_saints_in_paradise_lost():
    """
    Downloads, cleans, and searches the text of "Paradise Lost" for the names
    of historical saints.
    """
    print("Step 1: Downloading the text of 'Paradise Lost'...")
    # URL for the plain text version of Paradise Lost from Project Gutenberg
    url = "https://www.gutenberg.org/cache/epub/20/pg20.txt"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        full_text = response.text
        print("Download complete.")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. {e}")
        return

    print("Step 2: Cleaning the text (removing Gutenberg header/footer)...")
    try:
        # The main content is between these markers
        start_marker = "*** START OF THE PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        end_marker = "*** END OF THE PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        start_index = full_text.find(start_marker) + len(start_marker)
        end_index = full_text.find(end_marker)
        poem_text = full_text[start_index:end_index]
        print("Text cleaned successfully.")
    except Exception:
        print("Could not find start/end markers. Using full text, results may be inaccurate.")
        poem_text = full_text

    # Step 3: Define the list of historical saints to search for.
    # This list includes key apostles and major post-biblical saints.
    saints_to_find = [
        "Peter", "Paul", "Augustine", "Ambrose", "Jerome", "Gregory", "Francis",
        "Dominic", "Aquinas", "Benedict", "Athanasius", "Ignatius", "Polycarp"
    ]

    found_saints = set()

    print(f"Step 3: Searching for the names of {len(saints_to_find)} saints in the poem...")
    # Search for each saint in the cleaned text
    for saint in saints_to_find:
        # Use regex to find whole words only (\b), case-insensitively
        if re.search(r'\b' + re.escape(saint) + r'\b', poem_text, re.IGNORECASE):
            found_saints.add(saint)

    # Step 4: Display the results
    found_saints_list = sorted(list(found_saints))
    count = len(found_saints_list)

    print("\n--- Results ---")
    print("Historical saints mentioned by name in 'Paradise Lost':")
    if count == 0:
        print("None were found.")
    else:
        print(", ".join(found_saints_list))

    # Output the final equation as requested
    print("\nFinal Equation:")
    numbers_to_sum = ["1"] * count
    if not numbers_to_sum:
        # Equation for a zero result
        print("0 = 0")
    else:
        # Equation for a non-zero result
        print(f"{' + '.join(numbers_to_sum)} = {count}")


if __name__ == '__main__':
    count_saints_in_paradise_lost()