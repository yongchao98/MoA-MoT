def print_normal_cone_representation():
    """
    Prints the explicit representation of the normal cone T_F^°(x*).
    """
    
    print("The normal cone T_F^°(x^*) is the set of vectors s = (s1, s2, s3) in R^3.")
    print("The explicit representation is given by the following conditions on its components:")
    print("s1: any real number")
    print("s2: any real number")
    print("s3 >= 0")
    print("\nIn set-builder notation, this is:")
    print("T_F^°(x^*) = {s = (s1, s2, s3) in R^3 | s3 >= 0}")
    print("\nThis set represents the upper half-space in R^3, including the s1-s2 plane.")
    print("\nWe can express the condition s3 >= 0 using a dot product:")
    s_component = 1
    inequality_value = 0
    print(f"s^T * (0, 0, {s_component})^T >= {inequality_value}")

if __name__ == '__main__':
    print_normal_cone_representation()