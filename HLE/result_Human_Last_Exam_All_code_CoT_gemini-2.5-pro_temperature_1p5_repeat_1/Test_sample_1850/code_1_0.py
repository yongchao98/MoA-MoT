import requests
import re

def solve_milton_saints_puzzle():
    """
    Finds and counts the historical saints mentioned by name in Milton's Paradise Lost.
    """
    # Step 1 & 2: Acquire the text of Paradise Lost
    url = "https://www.gutenberg.org/files/20/20-0.txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the text file: {e}")
        return

    # Step 3: Define a dictionary of historical saints to search for.
    # We map alternate names to a canonical name to count each person only once.
    # This list includes key Old Testament figures venerated as saints.
    # New Testament saints like Peter or Paul are notably absent from Milton's text.
    saints_map = {
        "Adam": "Adam",
        "Eve": "Eve",
        "Abel": "Abel",
        "Enoch": "Enoch",
        "Noah": "Noah",
        "Abraham": "Abraham",
        "Abram": "Abraham",  # Alternate name for Abraham
        "Jacob": "Jacob",
        "Moses": "Moses",
        "Joshua": "Joshua",
        "David": "David",
        "Solomon": "Solomon",
    }

    found_saints = set()

    # Step 4: Search for each name and add the canonical name to a set.
    for search_name, canonical_name in saints_map.items():
        # Use regex to find whole words, ignoring case
        if re.search(r'\b' + re.escape(search_name) + r'\b', text, re.IGNORECASE):
            found_saints.add(canonical_name)

    # Step 5: Format and print the final output
    sorted_saints = sorted(list(found_saints))
    count = len(sorted_saints)

    print(f"There are {count} historical saints mentioned by name in Milton's Paradise Lost.")
    print("They are:")
    for saint in sorted_saints:
        print(f"- {saint}")

    # Create the final equation string
    equation_parts = ["1 (" + saint + ")" for saint in sorted_saints]
    equation_str = " + ".join(equation_parts)
    print("\nThe final equation is:")
    print(f"{equation_str} = {count}")


solve_milton_saints_puzzle()
<<<11>>>