import numpy as np

def solve_ququint_state():
    """
    Calculates the final state of a ququint after applying a given gate Q.
    """
    # 1. Define the 5x5 gate matrix Q from its action on the basis states.
    # Each column Q[:, j] represents the state Q|j>.
    Q = np.zeros((5, 5))
    sqrt2 = np.sqrt(2)

    # Q|0> = 1/sqrt(2) * (|1> + |2>)
    Q[1, 0] = 1 / sqrt2
    Q[2, 0] = 1 / sqrt2

    # Q|1> = 1/sqrt(2) * (|0> + |3>)
    Q[0, 1] = 1 / sqrt2
    Q[3, 1] = 1 / sqrt2

    # Q|2> = 1/sqrt(2) * (|1> + |4>)
    Q[1, 2] = 1 / sqrt2
    Q[4, 2] = 1 / sqrt2

    # Q|3> = 1/sqrt(2) * (|2> + |0>)
    Q[2, 3] = 1 / sqrt2
    Q[0, 3] = 1 / sqrt2

    # Q|4> = 1/sqrt(2) * (|3> + |2>)
    Q[3, 4] = 1 / sqrt2
    Q[2, 4] = 1 / sqrt2

    # 2. Define the initial state |ψ> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    psi = np.full(5, 1 / np.sqrt(5))

    # 3. Apply the gate Q to the state |ψ>
    psi_prime = Q @ psi

    # 4. The gate Q is not unitary, so the resulting state psi_prime is not normalized.
    # We will find the integer coefficients and then the correct normalization factor.
    # The common factor in psi_prime is 1/sqrt(10). Multiplying by sqrt(10) gives the integer coefficients.
    integer_coeffs = np.round(psi_prime * np.sqrt(10)).astype(int)

    # 5. The normalization constant is 1/sqrt(sum of squares of the integer coefficients).
    norm_squared = np.sum(integer_coeffs**2)

    # 6. Construct the final equation string for the normalized state.
    equation_parts = []
    for i, coeff in enumerate(integer_coeffs):
        if coeff != 0:
            # Handle the case where the coefficient is 1
            if coeff == 1:
                term = f"|{i}>"
            else:
                term = f"{coeff}|{i}>"
            equation_parts.append(term)

    final_equation = " + ".join(equation_parts)

    # Print the final state in the format: N * (c0|0> + c1|1> + ...)
    print("The final state of the system is:")
    print(f"|ψ'> = (1/sqrt({norm_squared})) * ({final_equation})")
    print("\nThis can also be written as:")
    
    # Print the equation with each coefficient explicitly shown
    full_equation_parts = []
    for i, coeff in enumerate(integer_coeffs):
        if coeff != 0:
            full_equation_parts.append(f"({coeff}/sqrt({norm_squared}))|{i}>")
    
    print("|ψ'> = " + " + ".join(full_equation_parts))


solve_ququint_state()
<<<|ψ'> = (1/sqrt(22)) * (2|0> + 2|1> + 3|2> + 2|3> + |4>)>>>