def solve_land_ownership():
    """
    This function determines and prints the final ownership of Lot A and Lot B
    based on the provided legal history.
    """

    # The final state of ownership is determined by the last legally recorded deeds.
    # Tommy's Deed was corrected and re-recorded for Lot B.
    # This makes Tommy's Estate the legal owner of Lot B.
    # The numbers are associated with the owner.
    tommy_estate_property = {
        "owner": "Tommy's Estate",
        "lot": "Lot B",
        "land_number": "1234567890",
        "trailer_number": "1234567890.1"
    }

    # Since Tommy's corrected deed is for Lot B, the intended swap
    # means James is the owner of Lot A, even though he never updated his paperwork.
    # His numbers are associated with him and therefore his property, Lot A.
    james_property = {
        "owner": "James",
        "lot": "Lot A",
        "land_number": "0987654321",
        "trailer_number": "0987654321.1"
    }

    # Print the final ownership details.
    print(f"{james_property['owner']} owns {james_property['lot']} with land number {james_property['land_number']} and trailer number {james_property['trailer_number']}.")
    print(f"{tommy_estate_property['owner']} owns {tommy_estate_property['lot']} with land number {tommy_estate_property['land_number']} and trailer number {tommy_estate_property['trailer_number']}.")

solve_land_ownership()

# The final answer is derived from the logic above.
# James owns Lot A with number 0987654321 and trailer 0987654321.1.
# Tommy's Estate owns Lot B with number 1234567890 and trailer 1234567890.1.