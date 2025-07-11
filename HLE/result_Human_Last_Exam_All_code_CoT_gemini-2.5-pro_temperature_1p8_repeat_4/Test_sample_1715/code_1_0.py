# The final answer needs to output each number in the final equation.
# This code will simulate the property ownership based on the legal events described.

# Represent the properties as dictionaries
tommy_property = {
    "owner": "Tommy",
    "lot": "A",
    "land_number": "1234567890",
    "trailer_number": "1234567890.1"
}

james_property = {
    "owner": "James",
    "lot": "B",
    "land_number": "0987654321",
    "trailer_number": "0987654321.1"
}

# In 2017, Tommy's deed was corrected and re-recorded.
# This is a legal change of ownership for the lot.
# Tommy's Estate now legally owns Lot B.
tommy_property["owner"] = "Tommy's Estate"
tommy_property["lot"] = "B" 

# James was advised to correct his deed to Lot A but never did.
# Legally, the intent was for him to have the other lot. 
# Although his paperwork is incorrect, the corrected deed for Tommy implies Lot A belongs to James.
james_property["lot"] = "A"

# Print the final ownership status based on the last recorded legal actions.
# We will present the information for Lot A first, then Lot B.

# Find who owns Lot A
owner_of_lot_a = None
if tommy_property["lot"] == "A":
    owner_of_lot_a = tommy_property
else:
    owner_of_lot_a = james_property

# Find who owns Lot B
owner_of_lot_b = None
if tommy_property["lot"] == "B":
    owner_of_lot_b = tommy_property
else:
    owner_of_lot_b = james_property

# Format and print the final answer
print("Based on the legally recorded deeds:")
print(f"Lot A is owned by {owner_of_lot_a['owner']} with Land Tax Map Number {owner_of_lot_a['land_number']} and Trailer number {owner_of_lot_a['trailer_number']}.")
print(f"Lot B is owned by {owner_of_lot_b['owner']} with Land Tax Map Number {owner_of_lot_b['land_number']} and Trailer number {owner_of_lot_b['trailer_number']}.")

final_answer = f"""
Based on the legally recorded deeds:
Lot A is owned by {owner_of_lot_a['owner']} with Land Tax Map Number {owner_of_lot_a['land_number']} and Trailer number {owner_of_lot_a['trailer_number']}.
Lot B is owned by {owner_of_lot_b['owner']} with Land Tax Map Number {owner_of_lot_b['land_number']} and Trailer number {owner_of_lot_b['trailer_number']}.
"""
print(f'<<<{final_answer}>>>')