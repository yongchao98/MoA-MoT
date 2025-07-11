import sys

def solve():
    """
    This program solves the bird identification task.
    
    1.  Analysis of the Target Image:
        The target Sociable Weaver displays a distinct scalloped feather pattern on its back. A key identifying feature is a small, dark irregularity or gap in the pattern on the upper right side of its back. The lighting, coloration, and overall arrangement of the feathers are also noted.
    
    2.  Comparison with Options:
        Each option (A-I) is compared against the target's features.
        - Options A, C, D, F, H, and I are eliminated due to differences in lighting, coloration, contrast, or visible pattern details.
        - Option B is a strong candidate with a similar pattern, but lacks the specific identifiable irregularity.
        - Option G perfectly matches the target image. It exhibits the same scalloped pattern, identical coloration, and most importantly, it shares the same small, dark irregularity on the upper right side of its back. The pose and lighting are also consistent.
    
    3.  Conclusion:
        Based on the matching of the unique feather pattern, especially the specific imperfection, image G is identified as the same bird as the target.
    """
    
    # The correct answer is determined by the detailed analysis above.
    correct_option = 'G'
    
    print(correct_option)

solve()