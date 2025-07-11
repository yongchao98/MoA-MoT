import datetime

# A list of notable French ships sunk by German U-boats (or their mines)
# between the start of the war and the armistice of 22 June 1940.
# GRT (Gross Register Tonnage) is used as a standard measure of a ship's size.
ships_sunk = [
    {'name': 'Pacifique', 'tonnage': 9474, 'sunk_date': datetime.date(1939, 9, 15)},
    {'name': 'Bretagne (cargo)', 'tonnage': 10108, 'sunk_date': datetime.date(1939, 10, 14)},
    {'name': 'Louisiane', 'tonnage': 6903, 'sunk_date': datetime.date(1939, 10, 17)},
    {'name': 'Flandre', 'tonnage': 8010, 'sunk_date': datetime.date(1940, 2, 13)},
    {'name': 'Brazza', 'tonnage': 10387, 'sunk_date': datetime.date(1940, 5, 28)},
    {'name': 'Champlain', 'tonnage': 28124, 'sunk_date': datetime.date(1940, 6, 17)},
]

# Initialize a variable to hold the largest ship found so far.
# Start with the first ship in the list.
largest_ship = ships_sunk[0]

print("Comparing the tonnage of French ships sunk by U-boats before the 1940 armistice:")
print("-" * 70)

# Iterate through the list of ships to find the one with the maximum tonnage.
for ship in ships_sunk:
    # Print each ship's details to show the comparison.
    print(f"Ship: {ship['name']:<20} Tonnage: {ship['tonnage']:>5} GRT")
    
    # Compare the current ship's tonnage with the largest one found so far.
    if ship['tonnage'] > largest_ship['tonnage']:
        largest_ship = ship

print("-" * 70)
print(f"The largest ship found is the {largest_ship['name']} with a tonnage of {largest_ship['tonnage']} GRT.")
print("\nThis passenger liner struck a mine laid by the German submarine U-65 off La Rochelle on June 17, 1940.")
print("While its GRT was 28,124, its full load displacement was approximately 30,500 tons, making it the largest French ship by displacement sunk by a U-boat before the armistice.")
