def find_common_elements():
    """
    This function identifies and prints the common literary elements between
    "Key West" by Koji Suzuki and "The Yellow Wallpaper" by Charlotte Perkins Gilman.
    """

    # Analysis based on a close reading of both texts.
    # The list below contains the Roman numerals for elements found to be in common.
    common_elements = [
        "I",   # Physical confinement
        "II",  # Progressive detachment from reality
        "III", # Indifference to others' emotions at the end
        "V",   # External locus of control
        "VI"   # Reference to a medical/biological conflict
    ]

    # The elements are already in numerical order.
    # Join the list with commas to create the final output string.
    final_answer = ", ".join(common_elements)
    
    print(final_answer)

find_common_elements()