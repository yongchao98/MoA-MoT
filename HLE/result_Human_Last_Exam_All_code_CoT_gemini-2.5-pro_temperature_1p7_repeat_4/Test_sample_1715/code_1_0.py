def solve_land_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the described deed correction.
    """
    # Initial state in 2016.
    # The numbers are tied to the person/property, not the Lot letter.
    tommy = {
        "owner": "Tommy",
        "tax_id": "1234567890",
        "trailer_id": "1234567890.1"
    }

    james = {
        "owner": "James",
        "tax_id": "0987654321",
        "trailer_id": "0987654321.1"
    }

    # In 2017, Tommy's Deed was corrected, legally giving his estate Lot B.
    # This implies the original intention was for Tommy to have Lot B and James to have Lot A.
    # A corrective deed supersedes the original incorrect one.

    final_lot_a_owner = james["owner"]
    final_lot_a_tax_id = james["tax_id"]
    final_lot_a_trailer_id = james["trailer_id"]

    final_lot_b_owner = "Tommy's Estate"
    final_lot_b_tax_id = tommy["tax_id"]
    final_lot_b_trailer_id = tommy["trailer_id"]
    
    # Print the final state of ownership for each lot.
    print(f"Lot A is owned by {final_lot_a_owner} with tax map number {final_lot_a_tax_id} and trailer number {final_lot_a_trailer_id}.")
    print(f"Lot B is owned by {final_lot_b_owner} with tax map number {final_lot_b_tax_id} and trailer number {final_lot_b_trailer_id}.")

solve_land_ownership()