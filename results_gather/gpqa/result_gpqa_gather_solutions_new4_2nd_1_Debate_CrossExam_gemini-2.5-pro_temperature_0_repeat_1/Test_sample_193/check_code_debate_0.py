import itertools
import collections
import sympy

def check_partition_function():
    """
    This function calculates the partition function for the given three-spin system
    and compares it to the provided answer.
    """
    # Define symbolic variables for the check
    J, beta = sympy.symbols('J beta')

    # 1. Enumerate all possible states
    spins = [-1, 1]
    # There are 2^3 = 8 states for (S1, S2, S3)
    states = list(itertools.product(spins, repeat=3))

    # 2. Calculate the energy for each state and find degeneracies
    energy_degeneracy = collections.Counter()
    for s1, s2, s3 in states:
        # Energy formula: E = -J[ S1S2 + S1S3 + S2S3 ]
        # The energy is expressed in units of J
        energy_coeff = -(s1 * s2 + s1 * s3 + s2 * s3)
        energy_degeneracy[energy_coeff] += 1

    # 3. Construct the partition function from the calculated degeneracies
    calculated_Z = 0
    for energy_coeff, degeneracy in energy_degeneracy.items():
        # The term is g * exp(-beta * E) = g * exp(-beta * energy_coeff * J)
        calculated_Z += degeneracy * sympy.exp(-beta * energy_coeff * J)

    # 4. Define the expression from the provided answer
    # The final answer is B, with the derived expression Z = 2e^(3Jβ) + 6e^(-Jβ)
    answer_expression_str = "2*exp(3*J*beta) + 6*exp(-J*beta)"
    try:
        # Use sympy's parser to create a symbolic expression from the string
        # This is safer than eval()
        answer_expression = sympy.sympify(answer_expression_str, locals={'J': J, 'beta': beta, 'exp': sympy.exp})
    except (sympy.SympifyError, SyntaxError):
        return "Failed to parse the expression from the provided answer."

    # 5. Check if the calculated expression matches the one in the answer
    # We use sympy.simplify to check for mathematical equivalence
    if sympy.simplify(calculated_Z - answer_expression) != 0:
        return (f"The derived expression in the answer is incorrect.\n"
                f"Calculated Z: {calculated_Z}\n"
                f"Answer's Z: {answer_expression}")

    # 6. Define the options from the question
    options = {
        'A': 2 * sympy.exp(-3*J*beta) + 6 * sympy.exp(J*beta),
        'B': 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta),
        'C': 2 * sympy.exp(2*J*beta) + 6 * sympy.exp(-2*J*beta),
        'D': 6 * sympy.exp(2*J*beta) + 2 * sympy.exp(-2*J*beta)
    }

    # The provided answer selected option B
    selected_option_letter = 'B'
    selected_option_expr = options.get(selected_option_letter)

    if selected_option_expr is None:
        return f"The selected option '{selected_option_letter}' is not a valid choice."

    # 7. Check if the selected option matches the correct calculation
    if sympy.simplify(calculated_Z - selected_option_expr) != 0:
        return (f"The final selected option '{selected_option_letter}' is incorrect.\n"
                f"The correct expression is Z = {calculated_Z}, but option {selected_option_letter} corresponds to Z = {selected_option_expr}.")

    return "Correct"

# Run the check
result = check_partition_function()
print(result)