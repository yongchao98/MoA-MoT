def solve_ownership():
    """
    This function tracks the ownership of Lot A and Lot B based on the provided legal events.
    """
    # Step 1: Define the initial ownership state in 2016.
    # Each person is represented by a dictionary containing their name, lot, and associated numbers.
    tommy = {
        "name": "Tommy",
        "lot": "A",
        "land_number": "1234567890",
        "trailer_number": "1234567890.1"
    }
    james = {
        "name": "James",
        "lot": "B",
        "land_number": "0987654321",
        "trailer_number": "0987654321.1"
    }

    # Step 2: Apply the changes that occurred in 2017.
    # Tommy's Deed was corrected and re-recorded from Lot A to Lot B.
    # His name is updated to reflect his estate's ownership.
    tommy["lot"] = "B"
    tommy["name"] = "Tommy's Estate"

    # Step 3: Determine the consequence for James.
    # Since Tommy's Estate now legally owns Lot B, James's claim on Lot B is superseded.
    # As advised by his lawyer, his ownership is now for Lot A.
    james["lot"] = "A"

    # Step 4: Identify who owns which lot for the final output.
    # We assign variables based on the final lot letter.
    if tommy["lot"] == "A":
        owner_A = tommy
        owner_B = james
    else:
        owner_A = james
        owner_B = tommy

    # Step 5: Print the final ownership details.
    print(f"Final Ownership of Lot A:")
    print(f"Owner: {owner_A['name']}")
    print(f"Lot: {owner_A['lot']}")
    print(f"Land Tax Map Number: {owner_A['land_number']}")
    print(f"Trailer Tax Map Number: {owner_A['trailer_number']}")
    print("\n" + "-"*30 + "\n")
    print(f"Final Ownership of Lot B:")
    print(f"Owner: {owner_B['name']}")
    print(f"Lot: {owner_B['lot']}")
    print(f"Land Tax Map Number: {owner_B['land_number']}")
    print(f"Trailer Tax Map Number: {owner_B['trailer_number']}")


solve_ownership()
<<<James owns Lot A with land number 0987654321 and trailer number 0987654321.1. Tommy's Estate owns Lot B with land number 1234567890 and trailer number 1234567890.1.>>>