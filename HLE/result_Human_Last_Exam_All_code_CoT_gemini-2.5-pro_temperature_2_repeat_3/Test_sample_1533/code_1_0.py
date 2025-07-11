def solve_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    Args:
        a (float): The length of side BC.
        b (float): The length of side AC.
        c (float): The length of side AB.
    
    Returns:
        float: The calculated ratio.
    """
    # Based on the geometric proof, the ratio BM / MI = (a + c) / b
    ratio = (a + c) / b
    return ratio

def main():
    """
    Main function to demonstrate the calculation for an example triangle.
    """
    # Example triangle side lengths (a scalene triangle)
    # a is the side opposite to vertex A (BC)
    # b is the side opposite to vertex B (AC)
    # c is the side opposite to vertex C (AB)
    a = 5.0
    b = 6.0
    c = 7.0

    # Calculate the ratio
    result = solve_ratio(a, b, c)

    # Print the derived formula and the final calculation
    print("The geometric analysis shows that the ratio can be expressed in terms of the side lengths a, b, and c.")
    print("The derived formula is: BM / MI = (a + c) / b")
    print("\nFor a triangle with the example side lengths:")
    print(f"a = {a} (side BC)")
    print(f"b = {b} (side AC)")
    print(f"c = {c} (side AB)")
    print("\nThe ratio is calculated as follows:")
    
    # Show each number in the final equation
    print(f"BM / MI = ({a} + {c}) / {b} = {result}")

if __name__ == "__main__":
    main()