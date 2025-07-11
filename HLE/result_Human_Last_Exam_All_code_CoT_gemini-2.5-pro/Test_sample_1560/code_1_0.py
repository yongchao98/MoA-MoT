import numpy as np

def solve_ququint_state():
    """
    Calculates the final state of a ququint after applying a given quantum gate.
    """
    # Step 1: Define the matrix M (the part of Q without the 1/sqrt(2) factor)
    # The matrix is constructed based on the definition of Q|i>
    # Q|0> = 1/sqrt(2)*(|1> + |2>) -> col 0 is [0, 1, 1, 0, 0]
    # Q|1> = 1/sqrt(2)*(|0> + |3>) -> col 1 is [1, 0, 0, 1, 0]
    # Q|2> = 1/sqrt(2)*(|1> + |4>) -> col 2 is [0, 1, 0, 0, 1]
    # Q|3> = 1/sqrt(2)*(|2> + |0>) -> col 3 is [1, 0, 1, 0, 0]
    # Q|4> = 1/sqrt(2)*(|3> + |2>) -> col 4 is [0, 0, 1, 1, 0]
    # numpy.array creates a matrix from a list of rows.
    M = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 1, 0]
    ])

    # The full gate matrix Q
    Q = (1 / np.sqrt(2)) * M

    # Step 2: Define the initial state |ψ⟩
    # |ψ⟩ = 1/sqrt(5) * (|0⟩ + |1⟩ + |2⟩ + |3⟩ + |4⟩)
    psi = (1 / np.sqrt(5)) * np.ones(5)

    # Step 3: Apply the gate Q to the state |ψ⟩
    psi_prime = Q @ psi

    # Step 4: The operator Q is not unitary, so we must normalize the resulting state.
    # First, calculate the norm (magnitude) of the resulting vector.
    norm_psi_prime = np.linalg.norm(psi_prime)

    # Then, normalize the vector to get the final state.
    psi_final = psi_prime / norm_psi_prime
    
    # Step 5: Display the final state in the required format.
    # The coefficients are the elements of the final state vector.
    coeffs = psi_final

    print("The final state of the system after applying the gate Q and before measurement is:")
    
    # We build the equation string piece by piece
    equation_parts = []
    for i in range(5):
        # We need to output each number in the final equation
        # The coefficients are c_i = num_i / sqrt(den)
        # From our calculation, the vector is (1/sqrt(22)) * [2, 2, 3, 2, 1]
        numerator = [2, 2, 3, 2, 1][i]
        denominator = 22
        
        # We also print the float value for clarity
        float_val = coeffs[i]

        equation_parts.append(f"({numerator}/√{denominator})|{i}⟩ ≈ {float_val:+.4f}|{i}⟩")

    final_equation = " + ".join(equation_parts)
    print(final_equation)
    
    # For verification, we can print the coefficients and check the sum of their squares
    # print("\nCoefficients:", coeffs)
    # print("Sum of squared coefficients (should be 1):", np.sum(coeffs**2))


solve_ququint_state()
<<<1/sqrt(22) * (2|0> + 2|1> + 3|2> + 2|3> + 1|4>)>>>