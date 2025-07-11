def solve_land_ownership():
    """
    This function determines and prints the final ownership of Lot A and Lot B
    based on the deed corrections described.
    """

    # The final state of ownership after the deed correction in 2017.
    # Tommy's deed was corrected to Lot B.
    # This implies a swap, leaving James with Lot A.
    # The tax/trailer numbers are tied to the person.

    final_ownership = {
        "Lot A": {
            "owner": "James",
            "number": "0987654321",
            "trailer_number": "0987654321.1"
        },
        "Lot B": {
            "owner": "Tommy's Estate",
            "number": "1234567890",
            "trailer_number": "1234567890.1"
        }
    }

    # Retrieve and print the details for Lot A
    lot_a_owner = final_ownership["Lot A"]["owner"]
    lot_a_number = final_ownership["Lot A"]["number"]
    lot_a_trailer = final_ownership["Lot A"]["trailer_number"]
    
    # Retrieve and print the details for Lot B
    lot_b_owner = final_ownership["Lot B"]["owner"]
    lot_b_number = final_ownership["Lot B"]["number"]
    lot_b_trailer = final_ownership["Lot B"]["trailer_number"]

    print("Based on the final recorded deeds:")
    print(f"{lot_a_owner} owns Lot A with the number {lot_a_number} and trailer number {lot_a_trailer}.")
    print(f"{lot_b_owner} owns Lot B with the number {lot_b_number} and trailer number {lot_b_trailer}.")

solve_land_ownership()
<<<Based on the final recorded deeds:
James owns Lot A with the number 0987654321 and trailer number 0987654321.1.
Tommy's Estate owns Lot B with the number 1234567890 and trailer number 1234567890.1.>>>