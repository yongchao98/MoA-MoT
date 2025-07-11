def calculate_riemann_components_kahler():
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.
    """
    try:
        m_str = input("Enter the complex dimension of the Kähler manifold (m): ")
        m = int(m_str)
        if m < 1:
            print("Error: The complex dimension must be a positive integer.")
            return

        print(f"\nFor a Kähler manifold of complex dimension m = {m}:")
        print("The real dimension n is 2 * m =", 2 * m)
        
        # The formula is (m * (m+1) / 2)^2
        # We will show each number in the equation
        
        m_plus_1 = m + 1
        numerator_1 = m * m_plus_1
        base = numerator_1 / 2
        
        # We need to handle the case where base is not an integer, but it always is
        # because either m or m+1 is even.
        base = int(base) 
        
        result = base ** 2

        print("\nThe formula for the number of independent components is (m * (m + 1) / 2)^2.")
        print("Let's substitute the value of m:")
        print(f"({m} * ({m} + 1) / 2)^2")
        print(f"= ({m} * {m_plus_1} / 2)^2")
        print(f"= ({numerator_1} / 2)^2")
        print(f"= ({base})^2")
        print(f"= {result}")
        
        print(f"\nResult: A Kähler manifold of complex dimension {m} has {result} independent Riemann tensor components.")

    except ValueError:
        print("Invalid input. Please enter an integer for the complex dimension.")

if __name__ == '__main__':
    calculate_riemann_components_kahler()