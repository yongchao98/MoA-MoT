def solve_beat_sheet_task():
    """
    This function identifies the indices of insect tribes whose immatures are
    unlikely to be collected using a beat-sheet method, based on their habitat.
    """
    
    # List of tribes provided
    tribes = [
        "1) Apis",
        "2) Melipotini",
        "3) Eupholini",
        "4) Acritini",
        "5) Oxyptilini",
        "6) Dictyophorini",
        "7) Acanthocerini"
    ]
    
    # Indices of tribes unlikely to be collected with a beat-sheet
    # 1: Apis larvae are in hives
    # 3: Eupholini larvae are typically root/stem borers
    # 4: Acritini larvae are predators in decaying matter
    # 7: Acanthocerini larvae are in dung balls underground
    unlikely_indices = [1, 3, 4, 7]
    
    # The problem asks for the indices in ascending order, separated by commas.
    # The list is already sorted, so we can directly format it for printing.
    result = ",".join(map(str, unlikely_indices))
    
    print(result)

solve_beat_sheet_task()