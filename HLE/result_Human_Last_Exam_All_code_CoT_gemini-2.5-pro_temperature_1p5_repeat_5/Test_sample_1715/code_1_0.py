def solve_ownership_puzzle():
    """
    This function tracks the ownership of Lot A and Lot B based on the provided story.
    It simulates the deed changes and prints the final result.
    """
    # Initial state in 2016 based on the first two deeds
    tommy = {
        "owner": "Tommy",
        "lot": "A",
        "land_id": "1234567890",
        "trailer_id": "1234567890.1"
    }

    james = {
        "owner": "James",
        "lot": "B",
        "land_id": "0987654321",
        "trailer_id": "0987654321.1"
    }

    # In 2017, Tommy's deed was corrected and re-recorded.
    # His ownership is now legally tied to Lot B.
    # This also means Lot A is now legally James's, even if his paperwork wasn't updated.
    
    # Update Tommy's record based on the re-recorded deed
    tommy["lot"] = "B"
    tommy["owner"] = "Tommy's Estate"
    
    # Update James's record to reflect the reciprocal change
    james["lot"] = "A"

    # Determine who owns which lot for the final output
    if tommy["lot"] == "A":
        owner_A = tommy
        owner_B = james
    else:
        owner_A = james
        owner_B = tommy

    # Print the final ownership details
    print("Based on the legally re-recorded deed:")
    print(f"Lot A is owned by: {owner_A['owner']}")
    print(f"Lot A Land Number: {owner_A['land_id']}")
    print(f"Lot A Trailer Number: {owner_A['trailer_id']}")
    print("-" * 20)
    print(f"Lot B is owned by: {owner_B['owner']}")
    print(f"Lot B Land Number: {owner_B['land_id']}")
    print(f"Lot B Trailer Number: {owner_B['trailer_id']}")

solve_ownership_puzzle()