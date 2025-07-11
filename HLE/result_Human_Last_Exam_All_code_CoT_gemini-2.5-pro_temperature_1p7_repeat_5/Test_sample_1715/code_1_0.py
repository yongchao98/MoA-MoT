# Initialize data for Tommy and James
tommy_estate = {
    "name": "Tommy's Estate",
    "land_id": "1234567890",
    "trailer_id": "1234567890.1",
    "lot_owned": "A"  # Initial state in 2016
}

james = {
    "name": "James",
    "land_id": "0987654321",
    "trailer_id": "0987654321.1",
    "lot_owned": "B"  # Initial state in 2016
}

# In 2017, Tommy's Deed for Lot A was found to be incorrect
# and was legally re-recorded for Lot B. This is the key legal change.
tommy_estate["lot_owned"] = "B"

# James was advised to correct his deed to Lot A, but he never did.
# Legally, the intended ownership for James is now Lot A,
# even though his paperwork is outdated. The correction on Tommy's deed
# effectively voids James's original deed for Lot B.
james_intended_lot = "A"


print("Current Ownership based on the latest recorded Deed:")
print("-" * 50)

# Output for Lot B
print(f"Lot {tommy_estate['lot_owned']} is legally owned by {tommy_estate['name']}.")
print(f"Associated Land Tax Map Number: {tommy_estate['land_id']}")
print(f"Associated Trailer Number: {tommy_estate['trailer_id']}")
print("-" * 50)

# Output for Lot A
print(f"Lot {james_intended_lot} is intended for {james['name']}, as his original Lot B deed was superseded by the correction.")
print("However, James's ownership is not legally recorded for Lot A because he never corrected his Deed.")
print(f"His original numbers are Land Tax Map Number: {james['land_id']} and Trailer Number: {james['trailer_id']}.")

final_answer = (
    f"Based on the re-recorded deed, {tommy_estate['name']} owns Lot {tommy_estate['lot_owned']} with numbers {tommy_estate['land_id']} and {tommy_estate['trailer_id']}. "
    f"James, who failed to correct his paperwork, physically occupies Lot {james_intended_lot} but his legal claim is unclear; his original deed for Lot {james['lot_owned']} is now in conflict with the corrected deed."
)
print("\n<<<")
print(final_answer)
print(">>>")
