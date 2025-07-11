import textdistance

def demonstrate_triangle_inequality():
    """
    This function demonstrates the triangle inequality for Levenshtein distance
    using the provided strings x, y, and z.
    The triangle inequality states: L(x, y) <= L(x, z) + L(z, y)
    """
    x = "algorithm"
    y = "logarithm"
    z = "altarithm"

    # Calculate the Levenshtein distances for each pair of strings
    l_xy = textdistance.levenshtein(x, y)
    l_xz = textdistance.levenshtein(x, z)
    l_zy = textdistance.levenshtein(z, y)

    # Print the result in the format of the inequality equation
    print("Verifying triangle inequality: L(x, y) <= L(x, z) + L(z, y)")
    print("For x = '{}', y = '{}', z = '{}'".format(x, y, z))
    
    # This fulfills the requirement to "output each number in the final equation"
    print("Result: {} <= {} + {}".format(l_xy, l_xz, l_zy))
    
    # Check if the inequality holds true and print a confirmation
    holds_true = l_xy <= l_xz + l_zy
    print("The inequality holds: {}".format(holds_true))

if __name__ == '__main__':
    demonstrate_triangle_inequality()
