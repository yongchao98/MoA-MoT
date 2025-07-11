def solve_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the legal changes described.
    """

    # Representing the individuals and their numbers
    tommy_owner = "Tommy's Estate"
    tommy_land_num = "1234567890"
    tommy_trailer_num = "1234567890.1"

    james_owner = "James"
    james_land_num = "0987654321"
    james_trailer_num = "0987654321.1"

    # In 2017, Tommy's Deed was legally corrected and re-recorded to Lot B.
    # This officially transfers ownership of Lot B to Tommy's Estate.
    final_lot_b_owner = tommy_owner
    final_lot_b_land_num = tommy_land_num
    final_lot_b_trailer_num = tommy_trailer_num

    # Since Tommy's Estate now legally owns Lot B, the remaining lot,
    # Lot A, belongs to James. James's failure to update his paperwork
    # does not change the legal status established by the re-recorded deed.
    final_lot_a_owner = james_owner
    final_lot_a_land_num = james_land_num
    final_lot_a_trailer_num = james_trailer_num

    # Print the final, legally correct ownership status
    print("Final Ownership Status:")
    print("-----------------------")
    print(f"Lot A is owned by {final_lot_a_owner}.")
    print(f"The associated land number is {final_lot_a_land_num} and the trailer number is {final_lot_a_trailer_num}.")
    print("\n")
    print(f"Lot B is owned by {final_lot_b_owner}.")
    print(f"The associated land number is {final_lot_b_land_num} and the trailer number is {final_lot_b_trailer_num}.")

solve_ownership()
<<<Tommy's Estate owns Lot B with land number 1234567890 and trailer number 1234567890.1. James owns Lot A with land number 0987654321 and trailer number 0987654321.1.>>>