def solve_ratio():
    """
    Calculates the ratio BM/MI for a triangle ABC given its side lengths a, b, c.
    
    The problem asks for the ratio BM/MI where:
    - I is the incenter of triangle ABC.
    - M is the intersection of the angle bisector BI and the circumcircle of ABC.
    
    The derived formula for this ratio is (a + c) / b.
    """
    
    print("Please provide the side lengths of the triangle ABC.")
    print("a is the length of the side opposite to vertex A (BC).")
    print("b is the length of the side opposite to vertex B (CA).")
    print("c is the length of the side opposite to vertex C (AB).")
    
    try:
        a = float(input("Enter side length a: "))
        b = float(input("Enter side length b: "))
        c = float(input("Enter side length c: "))
        
        # Triangle inequality check
        if a + b <= c or a + c <= b or b + c <= a:
            print("\nError: The given side lengths do not form a valid triangle.")
            return

        # Calculate the ratio using the derived formula
        ratio = (a + c) / b
        
        print("\nThe formula for the ratio BM/MI is (a + c) / b.")
        print("Substituting the given values:")
        # Final output showing the numbers in the equation
        print(f"({a} + {c}) / {b} = {ratio}")

    except ValueError:
        print("\nError: Please enter valid numerical values for the side lengths.")

if __name__ == '__main__':
    solve_ratio()
