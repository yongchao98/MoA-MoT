# Define the owners and their associated numbers
tommy_estate = {
    "name": "Tommy's Estate",
    "land_number": "1234567890",
    "trailer_number": "1234567890.1"
}

james = {
    "name": "James",
    "land_number": "0987654321",
    "trailer_number": "0987654321.1"
}

# The initial (incorrect) deed assignment
# Lot A -> Tommy
# Lot B -> James

# In 2017, Tommy's Deed was legally corrected and re-recorded.
# This makes the new assignment legally binding.
lot_b_owner = tommy_estate
lot_a_owner = james # By logical conclusion, as there are only two lots.

# Print the final ownership status based on the legal correction
# The question is "Who owns Lot A and Lot B with what number?"

print(f"Based on the legally re-recorded deed:")
print(f"Lot A is owned by {lot_a_owner['name']} with Tax Map Number {lot_a_owner['land_number']} and Trailer Number {lot_a_owner['trailer_number']}.")
print(f"Lot B is owned by {lot_b_owner['name']} with Tax Map Number {lot_b_owner['land_number']} and Trailer Number {lot_b_owner['trailer_number']}.")
print("\nNote: James's ownership of Lot A is based on the legal correction of Tommy's deed, even though James never updated his own paperwork.")
