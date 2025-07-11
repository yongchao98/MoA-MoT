import sys

def solve():
    """
    This function identifies the matching bird image.
    
    To re-identify the Sociable Weaver, I will analyze the features of the target image and compare them with the options provided.
    
    1.  **Analyze Target Image:** The target bird displays a distinct scalloped or scaled feather pattern on its back. The edges of these feathers are a bright, light grey, creating a high-contrast look. The overall coloration is greyish-brown. The bird is viewed from the back, angled slightly with its head turned to its right.
    
    2.  **Compare with Options:**
        -   **A, D, F:** These can be eliminated due to significant differences in lighting, angle, or overall coloration.
        -   **B, C, E, H, I:** These are closer, but subtle differences exist in the fineness of the pattern, the exact coloration, or the highlights.
        -   **G:** This image presents a bird with the same posture and angle as the target. The scalloped pattern on the back is a near-perfect match in terms of shape, size, and spacing. The lighting and coloration are also consistent with the target image. The way the light edges of the primary wing feathers are visible is identical.
    
    Therefore, image G shows the same bird as the target image.
    """
    # The letter corresponding to the correct image.
    correct_image_letter = 'G'
    print(correct_image_letter)

solve()