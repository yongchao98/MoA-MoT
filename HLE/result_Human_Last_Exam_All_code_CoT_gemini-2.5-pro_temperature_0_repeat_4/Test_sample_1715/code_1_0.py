def solve_land_ownership():
    """
    Analyzes the deed history to determine the current ownership of Lot A and Lot B.
    """

    # Define the numbers associated with each person
    tommy_land_num = "1234567890"
    tommy_trailer_num = "1234567890.1"

    james_land_num = "0987654321"
    james_trailer_num = "0987654321.1"

    # --- Event Analysis ---
    # In 2017, Tommy's Deed was corrected and re-recorded to Lot B.
    # This makes Tommy's Estate the legal owner of Lot B.
    final_lot_b_owner = "Tommy's Estate"
    final_lot_b_land_num = tommy_land_num
    final_lot_b_trailer_num = tommy_trailer_num

    # James was advised to correct his Deed to Lot A, but never did.
    # His recorded deed is still for Lot B. Therefore, based on the
    # provided documents, no one has a clear, recorded title to Lot A.
    final_lot_a_owner = "Disputed / Unclaimed"

    # --- Final Output ---
    print("Final Ownership Status Based on Recorded Deeds:")
    print("=" * 50)

    # Print Lot B details
    print("Lot B")
    print(f"  - Owner: {final_lot_b_owner}")
    print(f"  - Land Number: {final_lot_b_land_num}")
    print(f"  - Trailer Number: {final_lot_b_trailer_num}")
    print("  - Reason: Tommy's Deed was legally re-recorded to Lot B in 2017.")

    print("\n" + "=" * 50)

    # Print Lot A details
    print("Lot A")
    print(f"  - Owner: {final_lot_a_owner}")
    print("  - Reason: James never corrected his deed to Lot A.")
    print(f"  - Note: James's numbers ({james_land_num} and {james_trailer_num}) are tied to his uncorrected deed for Lot B.")

if __name__ == "__main__":
    solve_land_ownership()