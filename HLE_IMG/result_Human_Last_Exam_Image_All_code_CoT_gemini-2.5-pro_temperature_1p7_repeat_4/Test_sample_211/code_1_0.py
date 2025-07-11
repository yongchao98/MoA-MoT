import sys

def solve():
    """
    This function analyzes the provided images of Sociable weavers to find the match for the target image.

    1.  **Analyze Target Image:** The target bird displays a distinct scalloped pattern on its back feathers. Each feather has a dark center and a light, well-defined edge. The longer wing feathers have prominent light-colored stripes. The overall posture is a view from above, with the bird slightly angled. A key feature is a small, slightly darker or ruffled spot on the upper mantle, just below the neck area.

    2.  **Compare with Options:**
        *   **A, D, F, H, I:** These birds show different lighting, posture, or feather pattern details (e.g., sharpness, contrast, or arrangement of the 'scales') that do not perfectly align with the target.
        *   **B, C, E:** These are closer matches, but fine details in the pattern do not align perfectly when scrutinized.
        *   **G:** This image is a perfect match for the target. Every detail aligns:
            *   The posture is identical.
            *   The scalloped pattern on the back is an exact match, scale for scale.
            *   The prominent light-colored lines on the wings are in the same position and have the same thickness.
            *   The subtle imperfection (the small dark/ruffled spot) on the upper back is present in the exact same location.
            *   It is highly likely that the 'Target' image and image 'G' are crops from the same original photograph.

    Therefore, image G shows the same bird as the target image.
    """
    # The identified matching image is 'G'.
    answer = 'G'
    print(f"The image that shows the same bird as the target is: {answer}")

solve()