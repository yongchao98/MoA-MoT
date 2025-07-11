def solve_puzzle():
    """
    This function identifies the correct nyala based on visual features.
    
    1.  The target nyala has a unique cluster of white spots on its rear flank.
    2.  It also has a distinctive horn shape and stripe pattern.
    3.  Image A shows a younger animal with different horns and no spot cluster.
    4.  Image B is a good match, showing a similar spot cluster and horn shape.
    5.  Image C is of poor quality, making comparison difficult.
    6.  Image D is an exact match. The pattern of spots on the flank is identical to the target's. The horn shape is also a perfect match.
    7.  Image E has a different spot pattern and a marking on its leg.
    8.  Image F is obscured, preventing a clear comparison.
    9.  Comparing all options, Image D is the most definitive and clear match to the target nyala.
    """
    correct_option = 'D'
    print(f"The correct image is {correct_option}.")

solve_puzzle()