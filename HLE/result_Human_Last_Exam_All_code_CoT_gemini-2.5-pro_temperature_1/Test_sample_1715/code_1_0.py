def solve_land_ownership():
    """
    This function determines and prints the final ownership of Lot A and Lot B
    based on the described legal events.
    """

    # Initial state in 2016
    # tommy = {'owner': 'Tommy', 'lot': 'A', 'tax_id': '1234567890', 'trailer_id': '1234567890.1'}
    # james = {'owner': 'James', 'lot': 'B', 'tax_id': '0987654321', 'trailer_id': '0987654321.1'}

    # In 2017, Tommy's deed was corrected.
    # This is the decisive legal action.
    # Tommy's deed was re-recorded to Lot B.
    tommy_final = {'owner': "Tommy's Estate", 'lot': 'B', 'tax_id': '1234567890', 'trailer_id': '1234567890.1'}

    # Because Tommy's Estate now legally owns Lot B, James, by elimination and legal intent,
    # is the owner of Lot A, even though he failed to update his deed.
    james_final = {'owner': 'James', 'lot': 'A', 'tax_id': '0987654321', 'trailer_id': '0987654321.1'}

    print("Final Ownership Status:")
    print("-----------------------")
    
    # Print James's ownership details
    print(f"{james_final['owner']} owns Lot {james_final['lot']}.")
    print(f"The associated Tax Map Number is {james_final['tax_id']} and the trailer number is {james_final['trailer_id']}.")
    
    print("\n") # Add a newline for better separation

    # Print Tommy's Estate's ownership details
    print(f"{tommy_final['owner']} owns Lot {tommy_final['lot']}.")
    print(f"The associated Tax Map Number is {tommy_final['tax_id']} and the trailer number is {tommy_final['trailer_id']}.")


solve_land_ownership()