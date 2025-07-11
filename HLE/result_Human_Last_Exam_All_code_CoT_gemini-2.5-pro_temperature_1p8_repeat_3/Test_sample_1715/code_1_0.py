def solve_land_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the sequence of events described.
    """

    # Initial state in 2016
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
    # The owner's name is updated to reflect the estate.
    tommy["owner"] = "Tommy's Estate"
    tommy["lot"] = "B"

    # James was supposed to correct his deed to Lot A but never did.
    # However, since Tommy's Estate legally owns Lot B, Lot A is left for James
    # as per the original intent to swap the incorrectly filed lots.
    james["lot"] = "A"


    # Determine who owns which lot
    if tommy["lot"] == "A":
        owner_A = tommy
        owner_B = james
    else:
        owner_A = james
        owner_B = tommy

    print("Final Ownership Status:\n")

    print(f"{owner_A['owner']} owns Lot {owner_A['lot']}.")
    print(f"The associated Tax Map Number for the land is {owner_A['land_id']}.")
    print(f"The associated Tax Map Number for the trailer is {owner_A['trailer_id']}.\n")


    print(f"{owner_B['owner']} owns Lot {owner_B['lot']}.")
    print(f"The associated Tax Map Number for the land is {owner_B['land_id']}.")
    print(f"The associated Tax Map Number for the trailer is {owner_B['trailer_id']}.")

solve_land_ownership()
