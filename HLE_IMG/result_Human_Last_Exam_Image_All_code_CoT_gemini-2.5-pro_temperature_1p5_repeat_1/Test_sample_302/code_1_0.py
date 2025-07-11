def get_seal_characters():
    """
    This function identifies and prints the characters on the seal in simplified Chinese.
    The characters are read from top to bottom, right to left.
    The traditional characters are 長樂不老.
    The simplified characters are 长乐不老.
    """
    # The simplified Chinese characters corresponding to the seal script
    char1 = "长"  # Simplified of 長 (cháng)
    char2 = "乐"  # Simplified of 樂 (lè)
    char3 = "不"  # (bù) - same in traditional and simplified
    char4 = "老"  # (lǎo) - same in traditional and simplified

    # The characters are read from top-to-bottom, right-to-left.
    # So the order is char1, char2, char3, char4.
    seal_content = f"{char1}{char2}{char3}{char4}"
    
    print(f"The characters on the seal are: {seal_content}")
    print(f"In order: '{char1}', '{char2}', '{char3}', '{char4}'")

get_seal_characters()