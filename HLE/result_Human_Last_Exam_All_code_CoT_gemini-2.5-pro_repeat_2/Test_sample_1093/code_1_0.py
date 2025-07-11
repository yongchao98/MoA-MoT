import math

def solve_gauss_game():
    """
    Calculates the value that the average number of throws in Gauss's game
    converges to.

    This value is Euler's number, e, which can be expressed as the infinite
    sum of 1/k! for k from 0 to infinity.
    """
    # The number of terms to use for the approximation.
    # The series converges very quickly, so 18 is more than enough
    # for standard float precision.
    num_terms = 18
    
    approximation_e = 0.0
    
    print("The average of the recorded numbers converges to Euler's number, e.")
    print("The value of e is defined by the infinite sum: 1/0! + 1/1! + 1/2! + 1/3! + ...")
    print("\nHere is the calculation using the first few terms:")
    
    # Build and print the equation with its numerical values
    equation_parts = []
    for k in range(num_terms):
        term = 1.0 / math.factorial(k)
        approximation_e += term
        
        # We only display the first 8 terms to keep the output clean
        if k < 8:
            # Add the number 'term' to the equation string
            # The format {:.5f} shows the number with 5 decimal places
            equation_parts.append(f"{term:.5f}")

    # Join the parts with " + " and add "..." at the end
    equation_str = " + ".join(equation_parts) + " + ..."
    
    print(f"\ne ≈ {equation_str}")
    
    print(f"\nBy summing the first {num_terms} terms, we get a very precise approximation:")
    print(f"Average ≈ {approximation_e}")

# Run the solver
solve_gauss_game()