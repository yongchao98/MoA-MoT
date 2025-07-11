def solve_land_ownership():
    """
    This function determines the final ownership of Lot A and Lot B
    based on the sequence of recorded deeds.
    """
    
    # Initial state in 2016
    # The numbers are associated with the individuals and their assets (trailers/deeds).
    tommy = {
        "name": "Tommy's Estate",
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

    # In 2017, Tommy's deed for Lot A was found to be incorrect.
    # It was corrected and re-recorded, officially moving his ownership to Lot B.
    tommy["lot"] = "B"

    # The consequence of Tommy taking ownership of Lot B is that James's claim
    # on Lot B is superseded. Even though James never filed the correcting paperwork,
    # the intended swap leaves Lot A as his property.
    james["lot"] = "A"

    # Assign final ownership based on the "lot" key in each dictionary.
    if tommy["lot"] == "A":
        lot_a_owner = tommy
        lot_b_owner = james
    else:
        lot_a_owner = james
        lot_b_owner = tommy
        
    # Print the final, resolved ownership status.
    # The prompt asks for each number in the final equation.
    print(f"Based on the legally recorded deed correction:")
    print(f"{lot_a_owner['name']} owns Lot {lot_a_owner['lot']} with Tax Map Number {lot_a_owner['land_number']} and trailer number {lot_a_owner['trailer_number']}.")
    print(f"{lot_b_owner['name']} owns Lot {lot_b_owner['lot']} with Tax Map Number {lot_b_owner['land_number']} and trailer number {lot_b_owner['trailer_number']}.")

solve_land_ownership()
<<<James owns Lot A with Tax Map Number 0987654321 and trailer number 0987654321.1. Tommy's Estate owns Lot B with Tax Map Number 1234567890 and trailer number 1234567890.1.>>>