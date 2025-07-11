def solve_land_ownership():
    """
    Determines land ownership based on a series of recorded deeds.
    """

    # Initial state in 2016 based on the first deeds
    ownership_2016 = {
        "Tommy": {"Lot": "A", "Tax Number": "1234567890", "Trailer Number": "1234567890.1"},
        "James": {"Lot": "B", "Tax Number": "0987654321", "Trailer Number": "0987654321.1"}
    }

    # In 2017, Tommy's Deed was legally corrected and re-recorded.
    # This is the crucial event that changes the legal ownership.
    final_ownership = ownership_2016.copy()
    final_ownership["Tommy's Estate"] = final_ownership.pop("Tommy")
    final_ownership["Tommy's Estate"]["Lot"] = "B" # The deed was re-recorded to Lot B.

    # James never corrected his deed to claim Lot A.
    # His existing deed is for Lot B, which conflicts with Tommy's Estate's corrected deed.
    # No one has a filed and recorded deed for Lot A.

    tommy_estate_info = final_ownership["Tommy's Estate"]
    james_info = final_ownership["James"]

    print("Final Ownership Status based on recorded deeds:\n")

    # Ownership of Lot A
    print("Lot A:")
    print("No one holds a current, valid, recorded deed for Lot A.")
    print("Tommy's original deed for Lot A was nullified and re-recorded for Lot B.")
    print("James never filed a deed to claim Lot A.\n")


    # Ownership of Lot B
    print("Lot B:")
    print(f"Owned by: Tommy's Estate")
    print(f"This is based on the corrected and re-recorded Deed from 2017.")
    print(f"Associated Tax Map Number: {tommy_estate_info['Tax Number']}")
    print(f"Associated Trailer Number: {tommy_estate_info['Trailer Number']}\n")
    
    # Note on James's situation
    print("James's Situation:")
    print("James still holds an uncorrected 2016 deed and a lien for Lot B.")
    print(f"His original numbers are {james_info['Tax Number']} and {james_info['Trailer Number']}.")
    print("However, the 2017 re-recorded deed for Tommy's Estate likely supersedes his claim to Lot B, creating a legal conflict that needs resolution.")


solve_land_ownership()

# Final conclusion summarizing the most direct answer
answer = "Tommy's Estate owns Lot B (1234567890, 1234567890.1). Lot A is legally unclaimed as no one holds a valid recorded deed for it."
print(f"\n<<<Who owns Lot A and Lot B with what number?\n{answer}>>>")