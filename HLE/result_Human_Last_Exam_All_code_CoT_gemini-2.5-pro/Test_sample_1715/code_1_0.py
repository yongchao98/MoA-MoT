# This script determines the final ownership of Lot A and Lot B based on the provided story.

# The core of the problem is the deed correction in 2017.
# Initially in 2016:
# Tommy owned Lot A (Number: 1234567890, Trailer: 1234567890.1)
# James owned Lot B (Number: 0987654321, Trailer: 0987654321.1)

# In 2017, Tommy's Deed for Lot A was found to be incorrect and was re-recorded to Lot B.
# This action legally transfers ownership of Lot B to Tommy's Estate.
# The text states: "Now Tommy number is 1234567890 and trailer 1234567890.1 Lot B."

# Consequently, even though James never corrected his deed, the ownership of Lot A
# falls to him, as there cannot be two legal owners for Lot B.

# Final Ownership Status:
tommy_estate = {
    "owner": "Tommy's Estate",
    "lot": "Lot B",
    "land_number": "1234567890",
    "trailer_number": "1234567890.1"
}

james = {
    "owner": "James",
    "lot": "Lot A",
    "land_number": "0987654321",
    "trailer_number": "0987654321.1"
}

# Print the final results
print("Final Ownership Status:")
print("-----------------------")
print(f"{tommy_estate['owner']} owns {tommy_estate['lot']} with Land Tax Map Number {tommy_estate['land_number']} and Trailer Number {tommy_estate['trailer_number']}.")
print(f"{james['owner']} owns {james['lot']} with Land Tax Map Number {james['land_number']} and Trailer Number {james['trailer_number']}.")
