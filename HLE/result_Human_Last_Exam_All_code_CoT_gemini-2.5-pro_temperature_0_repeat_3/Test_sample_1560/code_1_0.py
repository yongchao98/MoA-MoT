import numpy as np

def solve_ququint_state():
    """
    This function calculates the final state of a ququint system after applying a given gate Q.
    """
    # Step 1: Define the gate Q as a matrix and the initial state |psi> as a vector.
    # The columns of the matrix Q are the results of Q acting on the basis states |0> through |4>.
    # The definitions are given as:
    # Q|0> = 1/sqrt(2) * (|1> + |2>)
    # Q|1> = 1/sqrt(2) * (|0> + |3>)
    # Q|2> = 1/sqrt(2) * (|1> + |4>)
    # Q|3> = 1/sqrt(2) * (|2> + |0>)
    # Q|4> = 1/sqrt(2) * (|3> + |2>)
    # We define the core matrix without the 1/sqrt(2) factor first.
    Q_core_matrix = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # The initial state is |psi> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>).
    # The core vector is [1, 1, 1, 1, 1].
    psi_core_vector = np.ones(5)

    # Step 2: Calculate the new state vector by applying the gate Q to |psi>.
    # The full operation is |psi'> = (1/sqrt(2) * Q_core) @ (1/sqrt(5) * psi_core).
    # This simplifies to |psi'> = (1/sqrt(10)) * (Q_core @ psi_core).
    # We calculate the resulting core vector first.
    psi_prime_core_vector = Q_core_matrix @ psi_core_vector

    # Step 3: Normalize the final state.
    # A physical state vector must have a norm of 1. We find the norm of our resulting
    # vector and divide by it. The norm of (c * v) is |c| * ||v||.
    # The squared norm of psi_prime_core_vector is the sum of the squares of its components.
    squared_norm = np.sum(psi_prime_core_vector**2)

    # Step 4: Construct and print the final equation for the normalized state.
    # The final state is the superposition before measurement. A measurement would collapse it
    # to a single basis state.
    print("The final state of the system, |psi_final>, before measurement is described by the equation:")

    equation_terms = []
    for i in range(len(psi_prime_core_vector)):
        numerator = int(psi_prime_core_vector[i])
        denominator_sqrt = int(squared_norm)
        # We construct each term of the form "(c)|i>"
        term = f"({numerator}/sqrt({denominator_sqrt}))|{i}>"
        equation_terms.append(term)

    # Join the terms with " + " to form the full equation string.
    final_equation = "|psi_final> = " + " + ".join(equation_terms)

    print(final_equation)

# Execute the function to solve the problem.
solve_ququint_state()