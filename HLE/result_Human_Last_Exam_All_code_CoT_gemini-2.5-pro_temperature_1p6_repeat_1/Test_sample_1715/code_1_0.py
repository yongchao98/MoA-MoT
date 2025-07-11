def solve_ownership_puzzle():
    """
    Analyzes and prints the final ownership status of Lot A and Lot B
    based on the sequence of recorded legal events.
    """

    # Final recorded status of the deeds after all events
    tommy_estate_claim = {
        "owner": "Tommy's Estate",
        "lot": "B",
        "land_number": "1234567890",
        "trailer_number": "1234567890.1",
        "deed_action": "Corrected and re-recorded in 2017"
    }

    james_claim = {
        "owner": "James",
        "lot": "B",
        "land_number": "0987654321",
        "trailer_number": "0987654321.1",
        "deed_action": "Original from 2016, never corrected to Lot A"
    }

    print("--- Final Ownership Status Based on Recorded Deeds ---")
    print("\n")

    # Conclusion for Lot A
    print("Lot A Ownership:")
    print("Based on the documents, no one has a valid, current claim to Lot A.")
    print("Tommy's Estate's claim was legally moved to Lot B, and James never corrected his deed to claim Lot A.")
    print("\n")

    # Conclusion for Lot B
    print("Lot B Ownership:")
    print("The title for Lot B is conflicted. The most recent corrective deed dictates ownership, but the old deed was never resolved.")
    print(f"1. {tommy_estate_claim['owner']} legally owns Lot B with Tax Map Number {tommy_estate_claim['land_number']} and Trailer Number {tommy_estate_claim['trailer_number']}. This is based on the 2017 corrective deed.")
    print(f"2. {james_claim['owner']} still holds an unresolved 2016 deed (and lien) for Lot B with Tax Map Number {james_claim['land_number']} and Trailer Number {james_claim['trailer_number']}, but it conflicts with the more recent corrective deed.")


if __name__ == "__main__":
    solve_ownership_puzzle()