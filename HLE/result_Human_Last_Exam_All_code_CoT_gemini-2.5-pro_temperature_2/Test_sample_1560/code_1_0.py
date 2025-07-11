import numpy as np

def solve_ququint_state():
    """
    This function calculates the final state of a ququint after applying a given quantum gate Q.
    """

    # 1. Define the matrix M (the part of Q without the 1/sqrt(2) factor).
    # The columns of M are derived from the definitions of Q|j>.
    # Q|0> = 1/sqrt(2)*(|1> + |2>) -> col 0 = [0, 1, 1, 0, 0]
    # Q|1> = 1/sqrt(2)*(|0> + |3>) -> col 1 = [1, 0, 0, 1, 0]
    # Q|2> = 1/sqrt(2)*(|1> + |4>) -> col 2 = [0, 1, 0, 0, 1]
    # Q|3> = 1/sqrt(2)*(|2> + |0>) -> col 3 = [1, 0, 1, 0, 0]
    # Q|4> = 1/sqrt(2)*(|3> + |2>) -> col 4 = [0, 0, 1, 1, 0]
    M = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # The gate Q matrix
    Q = (1 / np.sqrt(2)) * M

    # 2. Define the initial state vector |ψ⟩
    # |ψ⟩ = 1/sqrt(5) * (|0⟩ + |1⟩ + |2⟩ + |3⟩ + |4⟩)
    psi_vector = (1 / np.sqrt(5)) * np.ones(5)

    # 3. Calculate the final state vector |ψ'⟩ = Q|ψ⟩
    final_psi_vector = Q @ psi_vector

    # 4. Extract the coefficients for the final equation.
    # We can calculate the integer coefficients by multiplying M by a vector of ones.
    # The normalization factor will be (1/sqrt(2)) * (1/sqrt(5)) = 1/sqrt(10).
    unnormalized_coeffs = M @ np.ones(5)
    
    c0 = int(unnormalized_coeffs[0])
    c1 = int(unnormalized_coeffs[1])
    c2 = int(unnormalized_coeffs[2])
    c3 = int(unnormalized_coeffs[3])
    c4 = int(unnormalized_coeffs[4])
    
    # The denominator for the normalization constant C = 1/sqrt(norm_sq_den)
    norm_sq_den = 10
    
    # 5. Print the final state equation as requested.
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"The final state |ψ'⟩ after applying gate Q is:")
    print(f"|ψ'⟩ = (1/sqrt({norm_sq_den})) * ({c0}|0⟩ + {c1}|1⟩ + {c2}|2⟩ + {c3}|3⟩ + {c4}|4⟩)")


solve_ququint_state()
<<<|ψ'⟩ = (1/sqrt(10)) * (2|0⟩ + 2|1⟩ + 3|2⟩ + 2|3⟩ + 1|4⟩)>>>