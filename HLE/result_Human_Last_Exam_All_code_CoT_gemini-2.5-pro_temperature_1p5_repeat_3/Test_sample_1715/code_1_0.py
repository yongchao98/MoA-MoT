def solve_land_ownership():
    """
    Determines and prints the final ownership of Lot A and Lot B based on the provided legal history.
    """

    # Initial state in 2016
    tommy = {"name": "Tommy", "lot": "A", "tax_number": "1234567890", "trailer_number": "1234567890.1"}
    james = {"name": "James", "lot": "B", "tax_number": "0987654321", "trailer_number": "0987654321.1"}

    # In 2017, Tommy's Deed was legally corrected.
    # His ownership was re-recorded from Lot A to Lot B.
    # This legally supersedes the initial arrangement.
    tommy["lot"] = "B"

    # Because Tommy's ownership was corrected to Lot B, this implies the
    # initial assignment was a mistake and James's correct lot is Lot A.
    james["lot"] = "A"

    print("Based on the recorded deed corrections, here is the current ownership:")
    print("-" * 60)
    print(f"Owner: Tommy's Estate")
    print(f"Lot: {tommy['lot']}")
    print(f"Tax Map Number: {tommy['tax_number']}")
    print(f"Trailer Number: {tommy['trailer_number']}")
    print("-" * 60)
    print("Owner: James")
    print(f"Lot: {james['lot']} (pending his deed correction)")
    print(f"Tax Map Number: {james['tax_number']}")
    print(f"Trailer Number: {james['trailer_number']}")
    print("-" * 60)

    # Prepare the final answer string
    answer = (
        f"Tommy's Estate owns Lot {tommy['lot']} with number {tommy['tax_number']} and trailer {tommy['trailer_number']}. "
        f"James owns Lot {james['lot']} with number {james['tax_number']} and trailer {james['trailer_number']}."
    )
    
    print(f"\nFinal Answer: {answer}")


solve_land_ownership()
<<<Tommy's Estate owns Lot B with number 1234567890 and trailer 1234567890.1. James owns Lot A with number 0987654321 and trailer 0987654321.1.>>>