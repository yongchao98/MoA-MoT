import sys

def solve():
    """
    This function identifies the matching bird image.
    
    By closely examining the feather patterns, especially the unique scalloping on the back and the edging of the wing feathers, we can re-identify the individual bird.

    1.  **Target Bird Analysis:** The target bird has a distinct scalloped pattern on its back. There's a subtle irregularity in this pattern on the upper back, just below the head. The light-colored edges of the primary flight feathers are also a key feature.
    2.  **Comparison:**
        - Image A: Different lighting, pattern appears finer.
        - Image B: Coloration is much lighter.
        - Image D, E, F, G, H, I: Show variations in pattern, lighting, or angle that do not perfectly match the target's specific features.
        - Image C: The pose, lighting, and, most importantly, the specific scalloped pattern on the back, including the minor imperfection, are a near-perfect match to the target. The wing feather pattern also aligns.
    3.  **Conclusion:** Image C shows the same bird as the target image.
    """
    
    # The letter corresponding to the correct image
    correct_option = 'C'
    
    # Print the final answer
    print(correct_option)

solve()