def solve_square_counting():
    """
    Calculates the total number of squares on an n x n grid and
    displays the full summation formula.
    """
    # Set the size of the grid, n.
    n = 8

    # Initialize the total count and a list to store parts of the equation string.
    total_squares = 0
    equation_terms = []

    # The formula for the total number of squares is Sum_{m=1 to n} m^2 * (n-m+1).
    # Here, a = m and b = n - m + 1.
    # We loop from m = 1 to n to calculate each term of the sum.
    for m in range(1, n + 1):
        # 'a' and 'b' from the expression sum_{m=1 to n} a^2 * b
        a = m
        b = n - m + 1
        
        # Calculate the value of the current term.
        term_value = (a**2) * b
        
        # Add the term's value to the running total.
        total_squares += term_value
        
        # Store the string representation of the term, e.g., "1^2*8".
        equation_terms.append(f"{a}^2*{b}")

    # Construct the full equation string by joining the terms with " + ".
    full_equation = " + ".join(equation_terms)
    
    # Print the explanation and the final equation with all its numbers.
    print(f"For an {n}x{n} grid, the expressions are a=m and b=n-m+1.")
    print("The total number of squares is the sum over m from 1 to n:")
    print(f"Total = {full_equation}")
    
    # Print the final calculated result.
    print(f"\nTotal number of squares: {total_squares}")

solve_square_counting()