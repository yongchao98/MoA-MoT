def solve_nyala_identification():
    """
    This function analyzes and identifies the matching nyala based on key visual features.

    Analysis Steps:
    1.  Observe the Target Nyala: Key features are the lyre-shaped horns and a unique triangular pattern of three white spots on its rear flank.
    2.  Compare with Options:
        - A: Horns too small. Not a match.
        - B: Spot pattern on haunch does not match. Not a match.
        - C: Horns and coloration are different. Not a match.
        - D: Horn shape is a perfect match. The unique triangular spot pattern on the flank is identical to the target's. This is a definitive match.
        - E: Leg has a marking not on the target, and pose obscures key features. Not a likely match.
        - F: Key spot pattern is obscured by a bush. Cannot be confirmed.
    3.  Conclusion: Image D shows the same individual animal as the target image.
    """
    
    # The letter corresponding to the correct image.
    correct_answer = 'D'

    print(f"Based on the analysis of horn shape and the unique spot pattern, the correct image is: {correct_answer}")

# Execute the function to find and print the answer.
solve_nyala_identification()