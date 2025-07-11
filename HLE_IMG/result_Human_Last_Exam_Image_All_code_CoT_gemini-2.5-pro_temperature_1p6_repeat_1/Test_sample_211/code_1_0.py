def find_the_bird():
    """
    This function identifies the matching bird based on visual analysis.

    The identification process involves comparing the unique feather patterns,
    coloration, and posture of the target bird with the provided options.

    1.  **Target Analysis:** The target bird displays a distinct scalloped pattern
        on its back and wings. The head is turned sharply over its left shoulder.
    2.  **Comparison:** Each option is compared to the target.
        - Options A, B, D, E, F, G, H, and I show differences in feather patterns,
          lighting, coloration, or pose when compared to the target.
        - Option C shows a bird with an identical posture. The feather patterns on the
          back and wings, including the size, shape, and arrangement of the light
          edges, are a perfect match to the target.
    3.  **Conclusion:** The bird in image C is the same individual as the target.
    """
    # The correct option is 'C'.
    correct_option = 'C'
    print(f"The image that shows the same bird as the target is: {correct_option}")

find_the_bird()