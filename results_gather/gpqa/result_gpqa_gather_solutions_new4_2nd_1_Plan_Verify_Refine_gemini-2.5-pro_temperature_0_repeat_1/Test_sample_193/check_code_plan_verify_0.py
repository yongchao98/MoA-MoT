import itertools
from collections import Counter
import sympy

def check_partition_function():
    """
    Calculates the partition function for the given three-spin system and checks the correctness of the provided answer.

    The system has three spins S1, S2, S3, each taking values +1 or -1.
    The energy is E = -J[S1*S2 + S1*S3 + S2*S3].
    The partition function is Z = sum(exp(-beta * E)) over all states.

    The function will:
    1. Enumerate all 2^3 = 8 possible states.
    2. Calculate the energy for each state.
    3. Group states by energy to find the degeneracy of each energy level.
    4. Construct the symbolic partition function from the energies and degeneracies.
    5. Compare the derived partition function with the provided answer.
    """
    
    # Define symbolic variables for clarity in the final expression
    J, beta = sympy.symbols('J beta')

    # Each spin can be +1 or -1
    spins = [-1, 1]
    
    # Generate all 8 possible states (S1, S2, S3)
    all_states = list(itertools.product(spins, repeat=3))
    
    # Calculate the energy for each state and store them
    # We calculate E/J to get integer energy levels
    energy_levels_over_J = []
    for s1, s2, s3 in all_states:
        # Calculate the sum of products
        sum_of_products = s1*s2 + s1*s3 + s2*s3
        # Calculate the energy E = -J * (sum)
        # We store E/J for simplicity
        energy_over_J = -1 * sum_of_products
        energy_levels_over_J.append(energy_over_J)
        
    # Count the number of states for each energy level (degeneracy)
    degeneracy = Counter(energy_levels_over_J)
    
    # The expected degeneracies based on the correct calculation are:
    # E = -3J (E/J = -3) -> degeneracy = 2
    # E = +J  (E/J = +1) -> degeneracy = 6
    expected_degeneracy = {-3: 2, 1: 6}

    if degeneracy != expected_degeneracy:
        return f"Incorrect calculation of energy levels or degeneracies. Calculated degeneracies (for E/J) were {dict(degeneracy)}, but expected {expected_degeneracy}."

    # Construct the partition function Z from the calculated degeneracies
    calculated_Z = 0
    for energy_val_over_J, g in degeneracy.items():
        energy = energy_val_over_J * J
        calculated_Z += g * sympy.exp(-beta * energy)
        
    # The provided answer is B, which corresponds to Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    # Let's create the symbolic expression for option B
    correct_Z_expr = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)
    
    # The final answer from the LLM is 'B'
    final_answer_option = 'B'
    
    # Check if the calculated partition function matches the expression for option B
    # sympy.simplify helps ensure that expressions in different forms are compared correctly
    if sympy.simplify(calculated_Z - correct_Z_expr) == 0:
        # The calculation is correct and it matches option B.
        # Since the provided answer is B, the answer is correct.
        return "Correct"
    else:
        # This case should not be reached if the degeneracy check passes, but it's good practice.
        return f"The derived partition function {calculated_Z} does not match the expression for option {final_answer_option}, which is {correct_Z_expr}."

# Run the check
result = check_partition_function()
print(result)