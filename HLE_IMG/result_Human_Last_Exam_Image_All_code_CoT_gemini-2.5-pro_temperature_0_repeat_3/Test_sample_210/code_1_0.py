def solve_task():
    """
    This function identifies the correct nyala.
    The target nyala has a unique pattern of stripes and spots.
    By comparing these features, we can find the match.

    - Horn Shape: The target has long, lyre-shaped horns.
    - Stripe Pattern: Thin vertical white stripes.
    - Spot Pattern: A distinct cluster of white spots on the rear flank.

    Comparing these features with the options:
    - A: Horns are too small.
    - B: A good candidate, but the spot pattern in D is a more precise match.
    - C: Different build and coloration.
    - D: The horn shape, stripe pattern, and especially the unique spot pattern on the flank are a perfect match to the target.
    - E: Has a white bandage on its leg, which the target lacks.
    - F: Partially obscured, not a clear match.

    Therefore, D is the correct answer.
    """
    correct_image = 'D'
    print(f"The image showing the same nyala as the target is: {correct_image}")

solve_task()