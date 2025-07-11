import numpy as np

def solve_ququint_problem():
    """
    Solves the ququint state evolution problem.
    """
    # 1. Define the gate Q as a matrix.
    # The columns of the matrix are the results of applying Q to the basis vectors.
    # Q|0> = 1/sqrt(2) * (|1> + |2>) -> col 0 = [0, 1, 1, 0, 0]
    # Q|1> = 1/sqrt(2) * (|0> + |3>) -> col 1 = [1, 0, 0, 1, 0]
    # Q|2> = 1/sqrt(2) * (|1> + |4>) -> col 2 = [0, 1, 0, 0, 1]
    # Q|3> = 1/sqrt(2) * (|2> + |0>) -> col 3 = [1, 0, 1, 0, 0]
    # Q|4> = 1/sqrt(2) * (|3> + |2>) -> col 4 = [0, 0, 1, 1, 0]
    
    # We can represent Q by its unscaled integer part and the scaling factor.
    q_matrix_unscaled = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Define the initial state |ψ> as a vector.
    # |ψ> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    # We can represent its unscaled part as a vector of ones.
    psi_vector_unscaled = np.ones(5)

    # 3. Apply the unscaled gate to the unscaled state vector.
    # This gives the unscaled integer coefficients of the final state.
    # (Q_unscaled * 1/sqrt(2)) @ (psi_unscaled * 1/sqrt(5))
    # = (1/sqrt(10)) * (Q_unscaled @ psi_unscaled)
    final_vector_unscaled = q_matrix_unscaled @ psi_vector_unscaled

    # 4. Normalize the final state.
    # First, calculate the squared norm of the unscaled final vector.
    norm_squared = np.sum(final_vector_unscaled**2)
    
    # The coefficients of the normalized state are the unscaled coefficients
    # divided by the square root of the squared norm.
    
    # 5. Format and print the final equation.
    print("The final state of the system before measurement is described by the vector |ψ'>.")
    print("The calculation proceeds by applying gate Q to state |ψ> and normalizing the result.")
    
    terms = []
    for i in range(5):
        coefficient = int(final_vector_unscaled[i])
        if coefficient != 0:
            terms.append(f"({coefficient}/√{int(norm_squared)})|{i}⟩")

    final_equation = " + ".join(terms)
    
    print("\nThe final state is:")
    print(f"|ψ'⟩ = {final_equation}")
    
    # Return the final equation for the answer tag
    return final_equation

# Execute the function and capture the result for the final answer tag.
final_answer_equation = solve_ququint_problem()

# The final answer format required by the prompt
# The string replacement is to match the output style of the prompt more closely.
final_answer_equation = final_answer_equation.replace('√', 'sqrt(').replace('⟩','>') + ')'*(final_answer_equation.count('√'))
print(f"\n<<<|ψ'> = {final_answer_equation}>>>")
