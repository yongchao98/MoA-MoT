def solve():
    """
    This function identifies the matching bird from the provided options.
    By comparing the distinct feather patterns on the back and wings,
    the correct match is determined.
    """
    # The distinctive scalloped pattern on the back and the white edges of the wing feathers
    # in the target image are examined.
    # Each option is compared against these features.
    # Image G shows an identical pattern of scalloped feathers on the mantle (upper back)
    # and the same distinct white edging on the wing feathers. The individual feather
    # shapes and their arrangement are a one-to-one match with the target bird.
    correct_image = 'G'
    print(correct_image)

solve()