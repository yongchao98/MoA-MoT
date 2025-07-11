#
#
def solve_ownership():
    """
    This function tracks the ownership of Lot A and Lot B based on the provided story.
    It simulates the deed changes and prints the final ownership status.
    """

    # Initial state in 2016
    initial_ownership = {
        "Lot A": {"owner": "Tommy", "land_id": "1234567890", "trailer_id": "1234567890.1"},
        "Lot B": {"owner": "James", "land_id": "0987654321", "trailer_id": "0987654321.1"}
    }

    # In 2017, Tommy's Deed for Lot A was found to be incorrect.
    # It was re-recorded to Lot B. This is a legal correction.
    # This means Tommy's Estate now legally owns Lot B with Tommy's original numbers.
    # The original owner of Lot B (James) is consequently assigned Lot A.

    # Final state after the 2017 correction
    final_ownership = {
        "Lot B": {
            "owner": "Tommy's Estate",
            "land_id": initial_ownership["Lot A"]["land_id"],
            "trailer_id": initial_ownership["Lot A"]["trailer_id"]
        },
        "Lot A": {
            "owner": "James",
            "land_id": initial_ownership["Lot B"]["land_id"],
            "trailer_id": initial_ownership["Lot B"]["trailer_id"]
        }
    }

    # Print the final, legally-standing ownership details.
    lot_a_owner = final_ownership["Lot A"]["owner"]
    lot_a_land_id = final_ownership["Lot A"]["land_id"]
    lot_a_trailer_id = final_ownership["Lot A"]["trailer_id"]
    print(f"{lot_a_owner} owns Lot A with number {lot_a_land_id} and trailer {lot_a_trailer_id}.")

    lot_b_owner = final_ownership["Lot B"]["owner"]
    lot_b_land_id = final_ownership["Lot B"]["land_id"]
    lot_b_trailer_id = final_ownership["Lot B"]["trailer_id"]
    print(f"{lot_b_owner} owns Lot B with number {lot_b_land_id} and trailer {lot_b_trailer_id}.")


solve_ownership()
<<<James owns Lot A with number 0987654321 and trailer 0987654321.1. Tommy's Estate owns Lot B with number 1234567890 and trailer 1234567890.1.>>>