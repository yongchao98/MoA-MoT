def solve_land_ownership():
    """
    This function determines and prints the final ownership of Lot A and Lot B
    based on the legal changes described.
    """

    # The final ownership is determined by the last legally filed and recorded deed.
    # In 2017, Tommy's Deed was re-recorded to Lot B. This makes his estate the legal owner of Lot B.
    # Consequently, James is the owner of the remaining lot, Lot A.

    final_ownership = {
        "Tommy's Estate": {
            "Lot": "B",
            "Land_Number": "1234567890",
            "Trailer_Number": "1234567890.1"
        },
        "James": {
            "Lot": "A",
            "Land_Number": "0987654321",
            "Trailer_Number": "0987654321.1"
        }
    }

    # Print the final results clearly
    tommy_owner = "Tommy's Estate"
    tommy_info = final_ownership[tommy_owner]
    print(f"{tommy_owner} owns Lot {tommy_info['Lot']} with the number {tommy_info['Land_Number']} and trailer number {tommy_info['Trailer_Number']}.")

    james_owner = "James"
    james_info = final_ownership[james_owner]
    print(f"{james_owner} owns Lot {james_info['Lot']} with the number {james_info['Land_Number']} and trailer number {james_info['Trailer_Number']}.")

solve_land_ownership()

# The final answer summarized
final_answer = "Tommy's Estate owns Lot B (Number: 1234567890). James owns Lot A (Number: 0987654321)."
print(f"\n<<< {final_answer} >>>")