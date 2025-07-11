def solve_seal_characters():
    """
    This function identifies and prints the simplified Chinese characters
    from the provided seal image.
    The reading order is top-to-bottom, right-to-left.
    """
    # Characters identified from the seal, converted to simplified Chinese
    char_top_right = "马"
    char_bottom_right = "到"
    char_top_left = "功"
    char_bottom_left = "成"

    # The complete idiom
    idiom = char_top_right + char_bottom_right + char_top_left + char_bottom_left

    print("The characters on the seal in simplified Chinese are:")
    print(f"Top-Right: {char_top_right}")
    print(f"Bottom-Right: {char_bottom_right}")
    print(f"Top-Left: {char_top_left}")
    print(f"Bottom-Left: {char_bottom_left}")
    print(f"\nRead in the traditional order, they form the idiom: {idiom}")

solve_seal_characters()