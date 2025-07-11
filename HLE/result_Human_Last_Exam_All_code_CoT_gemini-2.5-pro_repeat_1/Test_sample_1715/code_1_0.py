def solve_land_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the provided legal deed history.
    """

    # Define the initial numbers associated with each person.
    tommy_land_number = '1234567890'
    tommy_trailer_number = '1234567890.1'

    james_land_number = '0987654321'
    james_trailer_number = '0987654321.1'

    # The 2017 re-recorded deed corrected Tommy's ownership from Lot A to Lot B.
    # This implies a swap, making James the owner of Lot A.

    # Final ownership of Lot A
    owner_A = "James"
    lot_A_land_number = james_land_number
    lot_A_trailer_number = james_trailer_number

    # Final ownership of Lot B
    owner_B = "Tommy's Estate"
    lot_B_land_number = tommy_land_number
    lot_B_trailer_number = tommy_trailer_number
    
    # Print the final status for each lot clearly.
    print(f"James owns Lot A with Tax Map Number {lot_A_land_number} and trailer number {lot_A_trailer_number}.")
    print(f"Tommy's Estate owns Lot B with Tax Map Number {lot_B_land_number} and trailer number {lot_B_trailer_number}.")

solve_land_ownership()