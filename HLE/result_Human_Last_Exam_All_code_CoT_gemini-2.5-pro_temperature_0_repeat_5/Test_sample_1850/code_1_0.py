import requests
import re

def count_saints_in_paradise_lost():
    """
    This function downloads Milton's Paradise Lost, searches for mentions of
    historical saints (interpreted as major Old and New Testament figures),
    and prints the names found and the total count.
    """
    # The URL for the plain text version of Paradise Lost
    url = "https://www.gutenberg.org/cache/epub/20/pg20.txt"

    # A list of historical/biblical figures revered as saints.
    saints_to_check = [
        "Abraham", "Andrew", "Bartholomew", "David", "Elijah", "Enoch",
        "Isaac", "Jacob", "James", "Job", "John", "Joseph", "Joshua",
        "Jude", "Lazarus", "Luke", "Mark", "Mary", "Matthew", "Matthias",
        "Moses", "Noah", "Paul", "Peter", "Philip", "Samson", "Simon",
        "Stephen", "Thaddeus", "Thomas"
    ]

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text from the URL. {e}")
        return

    # Clean the text to remove the Project Gutenberg header and footer
    poem_text = text
    try:
        start_marker = "*** START OF THIS PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        end_marker = "*** END OF THIS PROJECT GUTENBERG EBOOK PARADISE LOST ***"
        start_index = text.find(start_marker)
        end_index = text.find(end_marker)

        if start_index != -1 and end_index != -1:
            poem_text = text[start_index + len(start_marker):end_index]
    except Exception:
        # Fallback to using the whole text if parsing fails
        pass

    found_saints = []
    for saint in saints_to_check:
        # Use a regular expression to find the whole word, case-insensitive.
        pattern = r'\b' + re.escape(saint) + r'\b'
        if re.search(pattern, poem_text, re.IGNORECASE):
            found_saints.append(saint)

    count = len(found_saints)

    print(f"Found {count} historical saints (including Old Testament figures) mentioned by name:")
    
    if found_saints:
        # Create a simple "equation" for the final count
        equation_str = " + ".join(["1"] * count)
        for i, saint in enumerate(sorted(found_saints)):
            print(f"1 ({saint})")
        print(f"Final Equation: {equation_str} = {count}")
    else:
        print("No saints from the list were found.")
        print("Final Equation: 0 = 0")

# Execute the function
count_saints_in_paradise_lost()