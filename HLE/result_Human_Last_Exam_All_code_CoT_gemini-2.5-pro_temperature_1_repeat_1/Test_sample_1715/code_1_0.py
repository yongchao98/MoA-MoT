def solve_land_ownership():
    """
    This function tracks the ownership of Lot A and Lot B based on deed filings.
    """

    # Initial state in 2016 based on the first recorded deeds
    tommy = {
        "owner": "Tommy",
        "lot": "A",
        "tax_number": "1234567890",
        "trailer_number": "1234567890.1"
    }

    james = {
        "owner": "James",
        "lot": "B",
        "tax_number": "0987654321",
        "trailer_number": "0987654321.1"
    }

    # In 2017, Tommy's Deed was legally corrected and re-recorded.
    # This is a binding legal change that supersedes the original filing.
    # Tommy's deed was changed from Lot A to Lot B.
    tommy["lot"] = "B"
    tommy["owner"] = "Tommy's Estate"

    # The consequence of Tommy's deed being corrected to Lot B is that James's
    # ownership is now Lot A, even though he failed to update his paperwork.
    # The legal system corrected the initial swap.
    james["lot"] = "A"


    # Final Ownership Status
    print("Based on the legally re-recorded deed from 2017, here is the current ownership:")
    print("-" * 30)

    # Print Tommy's Estate's ownership
    final_tommy = tommy["owner"]
    final_tommy_lot = tommy["lot"]
    final_tommy_tax = tommy["tax_number"]
    print(f"{final_tommy} owns Lot {final_tommy_lot} with Tax Map Number {final_tommy_tax}.")

    # Print James's ownership
    final_james = james["owner"]
    final_james_lot = james["lot"]
    final_james_tax = james["tax_number"]
    print(f"{final_james} owns Lot {final_james_lot} with Tax Map Number {final_james_tax}.")
    print("-" * 30)
    print("\nNote: James should correct his deed paperwork for Lot B to Lot A to clear the title, as advised by his lawyer.")


solve_land_ownership()
<<<Tommy's Estate owns Lot B with number 1234567890. James owns Lot A with number 0987654321.>>>