# A list of notable French merchant ships sunk by U-boats before the June 22, 1940 armistice.
# Data includes the ship's name and its size in Gross Register Tonnage (GRT).
sunk_ships = [
    {'name': 'Bretagne', 'tonnage': 10125, 'date_sunk': '1939-10-14'},
    {'name': 'Champollion', 'tonnage': 12546, 'date_sunk': '1939-10-14'},
    {'name': 'Sheherazade', 'tonnage': 13467, 'date_sunk': '1940-06-21'},
    {'name': 'Flandre', 'tonnage': 8010, 'date_sunk': '1939-10-14'},
]

# Find the ship with the maximum tonnage from the list.
largest_ship = max(sunk_ships, key=lambda x: x['tonnage'])

# Print the result.
print("The largest French ship by displacement (tonnage) sunk by a U-boat before the 1940 armistice was the oil tanker SS Sheherazade.")
print("It was sunk on June 21, 1940, by the U-boat U-47.")
print("\nShip Details:")
print(f"Name: {largest_ship['name']}")
# The final part of the "equation" is the tonnage number itself.
print(f"Tonnage (GRT): {largest_ship['tonnage']}")