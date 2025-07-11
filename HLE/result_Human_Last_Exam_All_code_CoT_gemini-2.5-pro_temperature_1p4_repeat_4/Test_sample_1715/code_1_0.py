def solve_ownership_puzzle():
    """
    This function analyzes the property ownership based on the provided story
    and prints the final legal status of Lot A and Lot B.
    """

    # Represent the final, legally-recorded state of the properties after all events.
    # Tommy's Deed was re-recorded for Lot B in 2017. This is a definitive legal action.
    lot_b_owner = "Tommy's Estate"
    lot_b_deed_number = "1234567890"
    lot_b_trailer_number = "1234567890.1"

    # James's original deed was for Lot B. He was advised to correct it to Lot A but never did.
    # Therefore, he does not have a legally recorded deed for Lot A.
    # His numbers are tied to his original, now superseded, deed for Lot B.
    lot_a_intended_owner = "James"
    james_original_deed_number = "0987654321"
    james_trailer_number = "0987654321.1"

    print("Based on the legally filed and re-recorded documents, here is the final ownership status:")
    print("-" * 60)

    # Print the status of Lot B
    print("\nOwnership of Lot B:")
    print(f"The Deed for Tommy's property was corrected and re-recorded to Lot B in 2017.")
    print(f"This makes Tommy's Estate the legal owner of Lot B.")
    print(f"  Owner: {lot_b_owner}")
    print(f"  Tax Map Number: {lot_b_deed_number}")
    print(f"  Trailer Number: {lot_b_trailer_number}")

    print("-" * 60)

    # Print the status of Lot A
    print("\nOwnership of Lot A:")
    print(f"{lot_a_intended_owner} was advised to correct his deed to claim Lot A, but he never completed this action.")
    print(f"His only recorded deed is for Lot B, which is now legally owned by Tommy's Estate.")
    print(f"Therefore, {lot_a_intended_owner}'s ownership of Lot A is not legally established.")
    print(f"  Intended Owner: {lot_a_intended_owner}")
    print(f"  Associated Numbers (from invalid Lot B deed):")
    print(f"    Tax Map Number: {james_original_deed_number}")
    print(f"    Trailer Number: {james_trailer_number}")

if __name__ == '__main__':
    solve_ownership_puzzle()