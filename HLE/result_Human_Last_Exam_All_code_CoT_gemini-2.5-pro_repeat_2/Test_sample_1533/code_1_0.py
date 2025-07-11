def solve_ratio():
    """
    This function prints the derived ratio BM/MI in terms of the side lengths a, b, and c.
    """
    # The side lengths of the triangle ABC are denoted by a, b, and c.
    # a is the length of side BC
    # b is the length of side CA
    # c is the length of side AB
    
    # The derived expression for the ratio BM/MI is (a + c) / b.
    
    numerator_part_1 = 'a'
    numerator_operator = '+'
    numerator_part_2 = 'c'
    denominator = 'b'
    
    print("The ratio BM/MI can be expressed in terms of the side lengths a, b, and c as:")
    print(f"({numerator_part_1} {numerator_operator} {numerator_part_2}) / {denominator}")

solve_ratio()