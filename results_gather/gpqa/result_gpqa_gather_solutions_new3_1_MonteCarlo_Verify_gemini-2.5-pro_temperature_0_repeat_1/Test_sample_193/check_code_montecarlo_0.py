import sympy
import itertools
from collections import Counter

def check_partition_function():
    """
    This function calculates the partition function for the given 3-spin system
    and checks if the provided answer is correct.
    """
    # Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # The system has 3 spins, each can be +1 or -1.
    # Generate all 2^3 = 8 possible states.
    spins = [1, -1]
    all_states = list(itertools.product(spins, repeat=3))

    # Define the energy function
    def calculate_energy(state, J_var):
        s1, s2, s3 = state
        return -J_var * (s1 * s2 + s1 * s3 + s2 * s3)

    # Calculate the energy for each state
    energies = [calculate_energy(state, J) for state in all_states]

    # Count the degeneracy of each energy level
    energy_degeneracy = Counter(energies)

    # Construct the partition function Z from the calculated energies and degeneracies
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # The final answer provided by the LLM is C, which corresponds to:
    # C) Z = 2 * e^(3*J*β) + 6 * e^(-J*β)
    llm_answer_str = "C"
    
    # Define the expressions for all options symbolically
    options = {
        "A": 2 * sympy.exp(-3*J*beta) + 6 * sympy.exp(J*beta),
        "B": 6 * sympy.exp(2*J*beta) + 2 * sympy.exp(-2*J*beta),
        "C": 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta),
        "D": 2 * sympy.exp(2*J*beta) + 6 * sympy.exp(-2*J*beta)
    }

    # The expression corresponding to the LLM's chosen answer
    llm_expression = options.get(llm_answer_str)

    if llm_expression is None:
        return f"Invalid option '{llm_answer_str}' selected by the LLM."

    # Check if the calculated partition function matches the LLM's chosen expression
    # sympy.simplify(expr1 - expr2) == 0 is a robust way to check for symbolic equality
    if sympy.simplify(calculated_Z - llm_expression) == 0:
        # Further check the reasoning provided by the LLM
        # The LLM's reasoning correctly identifies two energy levels:
        # E1 = -3J with degeneracy g1 = 2
        # E2 = +J with degeneracy g2 = 6
        
        expected_degeneracy = {
            -3*J: 2,
            J: 6
        }
        
        if energy_degeneracy == expected_degeneracy:
            return "Correct"
        else:
            return (f"The final expression is correct, but the reasoning about energy levels is flawed. "
                    f"Calculated energy levels and degeneracies: {energy_degeneracy}. "
                    f"Expected: {expected_degeneracy}.")
    else:
        # Find which option, if any, the calculation matches
        correct_option = None
        for opt_key, opt_expr in options.items():
            if sympy.simplify(calculated_Z - opt_expr) == 0:
                correct_option = opt_key
                break
        
        return (f"Incorrect. The LLM chose option {llm_answer_str}, which corresponds to the expression {llm_expression}. "
                f"However, the correct partition function, calculated from first principles, is Z = {calculated_Z}. "
                f"This corresponds to option {correct_option}.")

# Run the check
result = check_partition_function()
print(result)