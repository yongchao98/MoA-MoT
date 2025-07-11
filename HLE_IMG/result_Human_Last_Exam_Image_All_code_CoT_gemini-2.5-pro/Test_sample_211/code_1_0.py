import sys

def solve():
    """
    This function identifies the matching bird image.
    
    To solve this, I will perform a visual comparison:
    1.  Observe the target image, focusing on the distinct scalloped feather pattern on the bird's back. Note the specific shapes, sizes, and the way the light reflects off them.
    2.  Pay close attention to any unique markers or imperfections. In the target image, there is a subtle irregularity in the pattern on the upper back, slightly to the left of the center.
    3.  Systematically compare this detailed pattern with each of the options (A-I).
    4.  Image A, C, D, E, F, H, and I can be eliminated due to significant differences in lighting, posture, or the overall flow of the feather pattern.
    5.  Image G is a very close contender, but the pattern on the upper back appears more uniform and lacks the specific subtle flaw seen in the target.
    6.  Image B shows a pattern, coloration, and lighting that are virtually identical to the target. Crucially, the subtle irregularity in the feather pattern on the upper back is also present in Image B, confirming it is the same individual.
    """
    # The letter corresponding to the correct image
    correct_image_letter = 'B'
    
    print(correct_image_letter)

solve()