import math

def solve_gauss_game(num_terms=10):
    """
    This function calculates the average number of throws for the described game.
    The theoretical answer is the mathematical constant e.
    This is calculated by summing the series 1/k! for k=0, 1, 2, ...

    The function will print each term in the sum as requested.
    """
    
    total_sum = 0.0
    equation_terms = []

    for k in range(num_terms):
        try:
            term = 1.0 / math.factorial(k)
            total_sum += term
            # Add the number from the term to the list for the final equation
            equation_terms.append(f"{term:.4f}")
        except OverflowError:
            print(f"Factorial of {k} is too large to compute, stopping here.")
            break

    # Format the equation string
    equation_string = " + ".join(equation_terms)

    print("The average number of throws converges to the constant 'e'.")
    print("This is calculated by the series: 1/0! + 1/1! + 1/2! + 1/3! + ...")
    print("\nThe equation showing the sum of the first {} terms is:".format(len(equation_terms)))
    # Print the final equation with each number
    print(f"{equation_string} + ... = {total_sum:.5f}")
    print(f"\nFor reference, the value of math.e is approximately {math.e:.5f}")


solve_gauss_game()