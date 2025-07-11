import pandas as pd

# Create a dictionary of Lady Butler's paintings and the years of the battles they depict.
# Format: { 'Painting Title': ('Battle/Event', Year) }
paintings_data = {
    'The Roll Call': ('Battle of Inkerman', 1854),
    'Quatre Bras': ('Battle of Quatre Bras', 1815),
    'Balaclava': ('Battle of Balaclava', 1854),
    'The Defence of Rorke\'s Drift': ('Battle of Rorke\'s Drift', 1879),
    'Scotland for Ever!': ('Battle of Waterloo', 1815),
    'Floreat Etona!': ('Battle of Laing\'s Nek', 1881),
    'The Camel Corps': ('Gordon Relief Expedition / Battle of Abu Klea', 1885),
    'The Remnants of an Army': ('1842 Retreat from Kabul', 1842)
}

# Find the painting depicting the latest battle.
# We use max() with a lambda function to find the item with the highest year in the dictionary values.
latest_painting, (latest_battle, latest_year) = max(paintings_data.items(), key=lambda item: item[1][1])

# Print the details for context
print("List of battles depicted by Lady Butler and their years:")
for painting, (battle, year) in paintings_data.items():
    print(f"- '{painting}' depicts the '{battle}' which occurred in {year}.")

print("\n--- Finding the latest battle ---")
# The "final equation" here is showing the components of the result.
print(f"The latest battle is from the painting: '{latest_painting}'.")
print(f"The battle is: '{latest_battle}'.")
print(f"The year of the battle is: {latest_year}.")
print(f"Therefore, the latest historical battle depicted is the {latest_battle} in the year {latest_year}.")
