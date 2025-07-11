def solve_land_ownership():
    """
    This function tracks the ownership of Lot A and Lot B based on recorded deeds.
    """
    # Initial state in 2016
    tommy = {
        'owner': "Tommy",
        'lot': "A",
        'land_num': "1234567890",
        'trailer_num': "1234567890.1"
    }
    james = {
        'owner': "James",
        'lot': "B",
        'land_num': "0987654321",
        'trailer_num': "0987654321.1"
    }

    # In 2017, Tommy's Deed was re-recorded, changing his ownership from Lot A to Lot B.
    # This action legally defines the new ownership status, regardless of James's inaction.
    tommy['lot'] = "B"
    tommy['owner'] = "Tommy's Estate"
    
    # As a consequence of the swap, James's property becomes Lot A.
    james['lot'] = "A"

    # Final Ownership Status
    # Sort the owners by Lot for clear presentation
    final_owners = sorted([tommy, james], key=lambda x: x['lot'])

    print("Based on the legally re-recorded deed in 2017, the current ownership is as follows:\n")

    for owner_data in final_owners:
      if owner_data['lot'] == 'A':
        james_owner = owner_data['owner']
        james_lot = owner_data['lot']
        james_land_num = owner_data['land_num']
        james_trailer_num = owner_data['trailer_num']
        print(f"{james_owner} owns {james_lot} with land number {james_land_num} and trailer number {james_trailer_num}.")
      elif owner_data['lot'] == 'B':
        tommy_owner = owner_data['owner']
        tommy_lot = owner_data['lot']
        tommy_land_num = owner_data['land_num']
        tommy_trailer_num = owner_data['trailer_num']
        print(f"{tommy_owner} owns {tommy_lot} with land number {tommy_land_num} and trailer number {tommy_trailer_num}.")


solve_land_ownership()
<<<James owns Lot A with land number 0987654321 and trailer number 0987654321.1. Tommy's Estate owns Lot B with land number 1234567890 and trailer number 1234567890.1.>>>