import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension for the described fractal piano.
    """
    # N is the number of self-similar copies. In one octave, there are 5 black keys,
    # and each is replaced by a smaller version of the whole keyboard.
    N = 5

    # r is the scaling factor.
    # The original keyboard is 3 units wide.
    # There are 7 white keys, so each is 3/7 units wide.
    # A black key is half as wide: (1/2) * (3/7) = 3/14 units.
    # To fit the 3-unit-wide keyboard into the 3/14-unit-wide black key space,
    # the scaling factor 'r' is (3/14) / 3 = 1/14.
    r = 1/14
    
    # The dimension D is log(N) / log(1/r)
    one_over_r = 1/r
    dimension = math.log(N) / math.log(one_over_r)

    print("The Minkowski–Bouligand dimension (D) for a self-similar set is calculated using the formula:")
    print("D = log(N) / log(1/r)\n")
    print(f"1. N, the number of non-overlapping copies, is the number of black keys in an octave. So, N = {N}.")
    print(f"2. r, the scaling factor, is the ratio of the new keyboard's width to the old one's. This is (3/14) / 3 = 1/14. So, 1/r = {one_over_r:.0f}.")
    print("\nPlugging the values into the formula:")
    print(f"D = log({N}) / log({one_over_r:.0f})")
    print(f"D = {dimension}")

calculate_fractal_dimension()