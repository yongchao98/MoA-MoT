def guess_the_country():
    """
    Analyzes a small part of a flag to identify the country.
    """
    # 1. The image shows three horizontal stripes in the top-left corner.
    colors = ["black", "white", "red"]

    # 2. The order of colors from top to bottom is Black, White, Red.
    print("The flag in the image has the following horizontal stripes from top to bottom:")
    print("- " + colors[0])
    print("- " + colors[1])
    print("- " + colors[2])
    print("\n")

    # 3. No current UN member state has a black-white-red horizontal tricolor.
    # 4. However, the flag of Yemen is a red-white-black horizontal tricolor.
    yemen_flag_colors = ["red", "white", "black"]
    print("The flag of Yemen has the following horizontal stripes from top to bottom:")
    print("- " + yemen_flag_colors[0])
    print("- " + yemen_flag_colors[1])
    print("- " + yemen_flag_colors[2])
    print("\n")

    # 5. The flag in the image is an exact inversion of the flag of Yemen.
    # This is a common error in flag quizzes.
    print("Observation: The flag shown is a vertical inversion of the flag of Yemen.")
    print("Conclusion: Despite the error in the image, the intended country is most likely Yemen.")

guess_the_country()