def solve_ownership_puzzle():
    """
    This function tracks the ownership of Lot A and Lot B based on the provided story
    and prints the final ownership status according to the legally recorded deeds.
    """

    # Step 1: Define the initial state of ownership in 2016.
    # The tax map numbers are associated with the deeds of the respective owners.
    initial_ownership = {
        "Tommy": {
            "lot": "A",
            "land_number": "1234567890",
            "trailer_number": "1234567890.1"
        },
        "James": {
            "lot": "B",
            "land_number": "0987654321",
            "trailer_number": "0987654321.1"
        }
    }

    # Step 2: Simulate the legal correction in 2017.
    # Tommy's Deed for Lot A was incorrect and was re-recorded to Lot B.
    # This means Tommy's Estate now legally owns Lot B.
    # The numbers follow the deed, so Tommy's numbers now apply to Lot B.
    
    final_ownership = {}

    # Lot B's ownership is transferred to Tommy's Estate via the corrected deed.
    final_ownership["B"] = {
        "owner": "Tommy's Estate",
        "land_number": initial_ownership["Tommy"]["land_number"],
        "trailer_number": initial_ownership["Tommy"]["trailer_number"]
    }

    # James never corrected his deed to claim Lot A. Legally, his ownership of Lot A is not recorded.
    # His original numbers were intended for the property he would own.
    final_ownership["A"] = {
        "owner": "James (Note: Deed was not corrected, so ownership is not legally recorded)",
        "land_number": initial_ownership["James"]["land_number"],
        "trailer_number": initial_ownership["James"]["trailer_number"]
    }
    
    # Step 3: Print the final conclusion clearly.
    
    print("--- Final Ownership Status Based on Recorded Deeds ---")
    
    # Details for Lot B
    lot_b = final_ownership["B"]
    print("\nLot B:")
    print(f"  Owner: {lot_b['owner']}")
    print(f"  Land Tax Map Number: {lot_b['land_number']}")
    print(f"  Trailer Tax Map Number: {lot_b['trailer_number']}")

    # Details for Lot A
    lot_a = final_ownership["A"]
    print("\nLot A:")
    print(f"  Owner Status: {lot_a['owner']}")
    print(f"  Associated Land Tax Map Number: {lot_a['land_number']}")
    print(f"  Associated Trailer Tax Map Number: {lot_a['trailer_number']}")


solve_ownership_puzzle()