def solve_seal_characters():
    """
    This function identifies the characters on the seal, explains the reading order,
    and prints the result in simplified Chinese.
    """
    # The characters are identified as 景, 園, 不, 翁 in traditional script.
    # The standard reading order for a 2x2 seal is top-right, bottom-right, top-left, bottom-left.
    
    # Traditional characters based on identification
    char_tr_tr = "景"  # Top-Right
    char_br_tr = "園"  # Bottom-Right
    char_tl_tr = "不"  # Top-Left
    char_bl_tr = "翁"  # Bottom-Left
    
    # Simplified Chinese conversion
    char_tr_sm = "景"
    char_br_sm = "园"
    char_tl_sm = "不"
    char_bl_sm = "翁"
    
    # The phrase in simplified Chinese
    result = char_tr_sm + char_br_sm + char_tl_sm + char_bl_sm
    
    print(f"The four characters on the seal are read in the order of top-right, bottom-right, top-left, and bottom-left.")
    print(f"The traditional characters are {char_tr_tr}{char_br_tr}{char_tl_tr}{char_bl_tr}.")
    print(f"The simplified Chinese characters are: {result}")

solve_seal_characters()