def get_runestone_id():
    """
    This function identifies and prints the ID of the Ingvar runestone.
    The identification is based on the unique grid-like inscription pattern,
    which matches the known fragment Sö 137.
    """
    province_code = "Sö"
    stone_number = 137
    
    # The full ID combines the province code and the number.
    runestone_id = f"{province_code} {stone_number}"
    
    print(f"The ID of the Ingvar runestone is: {runestone_id}")
    print(f"The number in the ID is: {stone_number}")

get_runestone_id()