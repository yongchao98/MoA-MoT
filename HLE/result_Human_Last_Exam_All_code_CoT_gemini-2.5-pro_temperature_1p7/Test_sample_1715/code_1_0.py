def solve_land_ownership():
    """
    This function tracks the ownership of Lot A and Lot B
    through the events described and prints the final status.
    """
    # Initial ownership in 2016
    ownership = {
        "1234567890": {"owner": "Tommy", "lot": "A", "trailer_id": "1234567890.1"},
        "0987654321": {"owner": "James", "lot": "B", "trailer_id": "0987654321.1"}
    }

    # In 2017, Tommy's Deed was corrected and re-recorded.
    # The property identified by his number is now legally Lot B.
    # Tommy's Estate inherits this corrected ownership.
    ownership["1234567890"]["lot"] = "B"
    ownership["1234567890"]["owner"] = "Tommy's Estate"

    # James was supposed to correct his deed to Lot A but didn't.
    # His paperwork is incorrect, but his property number is associated with the
    # parcel that is physically and legally designated as Lot A, following the correction.
    # Therefore, James owns Lot A.
    ownership["0987654321"]["lot"] = "A"

    # Extracting final information for clarity
    tommy_estate_property = ownership["1234567890"]
    james_property = ownership["0987654321"]

    # Print the final resolved ownership details
    print("Based on the legally corrected and re-recorded deed:")
    print("-" * 30)

    # Details for Lot A
    print(f"The owner of Lot {james_property['lot']} is {james_property['owner']}.")
    print(f"The Tax Map Number for the land is {list(ownership.keys())[1]}.")
    print(f"The Tax Map Number for the trailer is {james_property['trailer_id']}.")

    print("-" * 30)

    # Details for Lot B
    print(f"The owner of Lot {tommy_estate_property['lot']} is {tommy_estate_property['owner']}.")
    print(f"The Tax Map Number for the land is {list(ownership.keys())[0]}.")
    print(f"The Tax Map Number for the trailer is {tommy_estate_property['trailer_id']}.")
    print("-" * 30)

solve_land_ownership()
<<<Tommy's Estate owns Lot B with numbers 1234567890 and 1234567890.1. James owns Lot A with numbers 0987654321 and 0987654321.1.>>>