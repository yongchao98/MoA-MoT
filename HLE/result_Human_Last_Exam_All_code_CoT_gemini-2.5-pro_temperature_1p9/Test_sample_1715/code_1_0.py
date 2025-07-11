def solve_land_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the described sequence of events.
    """

    # Person data with their unique tax map numbers
    people = {
        "Tommy": {
            "name": "Tommy's Estate",
            "land_number": "1234567890",
            "trailer_number": "1234567890.1"
        },
        "James": {
            "name": "James",
            "land_number": "0987654321",
            "trailer_number": "0987654321.1"
        }
    }

    # In 2017, Tommy's Deed was found to be incorrect.
    # It was re-filed and recorded to change his ownership from Lot A to Lot B.
    # This corrected deed is the most current and legally binding document for his property.
    
    tommy_final_lot = "B"
    james_final_lot = "A" # By logical deduction in a two-lot split

    # Assign final ownership based on the corrected deed.
    final_ownership = {
        "Lot A": people["James"],
        "Lot B": people["Tommy"]
    }
    
    owner_lot_b = final_ownership["Lot B"]
    owner_lot_a = final_ownership["Lot A"]
    
    print("Final Ownership Status:")
    # Print the ownership for Lot B
    print(
        f"{owner_lot_b['name']} owns Lot B with Tax Map Number {owner_lot_b['land_number']} and trailer number {owner_lot_b['trailer_number']}."
    )
    
    # Print the ownership for Lot A
    print(
        f"{owner_lot_a['name']} owns Lot A with Tax Map Number {owner_lot_a['land_number']} and trailer number {owner_lot_a['trailer_number']}."
    )
    
    print("\nNote: The correction of Tommy's deed to Lot B in 2017 would legally take precedence over James's older, uncorrected 2016 deed for Lot B.")


solve_land_ownership()