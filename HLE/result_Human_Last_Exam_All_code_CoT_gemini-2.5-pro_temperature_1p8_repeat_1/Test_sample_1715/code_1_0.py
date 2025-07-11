# This script determines the final ownership of Lot A and Lot B
# based on the provided sequence of events.

# Initial state in 2016 after the first filing.
tommy_estate = {
    "owner": "Tommy's Estate",
    "lot": "Lot A",  # Initially recorded as Lot A
    "land_number": "1234567890",
    "trailer_number": "1234567890.1"
}

james = {
    "owner": "James",
    "lot": "Lot B", # Initially recorded as Lot B
    "land_number": "0987654321",
    "trailer_number": "0987654321.1"
}

# In 2017, a correction was made to Tommy's Deed because it was filed incorrectly.
# The re-recorded deed correctly assigns Lot B to Tommy's Estate.
# This is the legally binding document for Tommy's property.
tommy_estate["lot"] = "Lot B"

# Because Tommy's Estate's ownership of Lot B is now legally established,
# James's property is, by necessity, Lot A. Even though his deed was not
# updated, the correction on the other deed clarifies his true ownership.
james["lot"] = "Lot A"

# Print the final, correct ownership status.
print(f"{tommy_estate['owner']} owns {tommy_estate['lot']} with Tax Map Number {tommy_estate['land_number']} and Trailer Number {tommy_estate['trailer_number']}.")
print(f"{james['owner']} owns {james['lot']} with Tax Map Number {james['land_number']} and Trailer Number {james['trailer_number']}.")
