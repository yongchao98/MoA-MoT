def solve_land_ownership():
    """
    This function determines the current ownership of Lot A and Lot B
    based on the described legal changes.
    """

    # Initial numbers assigned to each person
    tommy_land_num = "1234567890"
    tommy_trailer_num = "1234567890.1"

    james_land_num = "0987654321"
    james_trailer_num = "0987654321.1"

    # In 2017, Tommy's Deed was corrected and re-recorded to Lot B.
    # This is the last and legally binding action.
    tommy_estate_owner = {
        "owner": "Tommy's Estate",
        "lot": "B",
        "land_num": tommy_land_num,
        "trailer_num": tommy_trailer_num
    }

    # Since the corrected deed gives Lot B to Tommy's Estate,
    # Lot A is left for James, as was advised by the lawyer.
    # James's failure to update his deed doesn't change the fact that
    # the deed for Lot B was legally transferred to Tommy's Estate.
    james_owner = {
        "owner": "James",
        "lot": "A",
        "land_num": james_land_num,
        "trailer_num": james_trailer_num
    }
    
    # Final ownership status
    print(f"{james_owner['owner']} owns Lot {james_owner['lot']} with Tax Map Number {james_owner['land_num']} and trailer number {james_owner['trailer_num']}.")
    print(f"{tommy_estate_owner['owner']} owns Lot {tommy_estate_owner['lot']} with Tax Map Number {tommy_estate_owner['land_num']} and trailer number {tommy_estate_owner['trailer_num']}.")

solve_land_ownership()