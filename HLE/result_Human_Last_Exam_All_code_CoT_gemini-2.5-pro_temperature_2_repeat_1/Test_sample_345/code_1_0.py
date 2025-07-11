import sys
import io

# Set stdout to utf-8 to handle special characters if any
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Based on art history and museum collection data, the 1927 tempera painting
# "The Radionist" ("Der Radionist") by Kurt Günther was acquired in 1967.
# The acquiring institution is a German museum.

museum_name = "Lindenau-Museum Altenburg"
painting_title = "The Radionist"
artist_name = "Kurt Günther"
year_of_acquisition = 1967

# Print the final answer
print(f"The {museum_name} acquired the painting '{painting_title}' by {artist_name} in the year {year_of_acquisition}.")
