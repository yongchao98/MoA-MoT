import sympy
from itertools import product
from collections import Counter

def check_partition_function_correctness():
    """
    This function programmatically calculates the partition function for a system
    of three spins and verifies it against a candidate answer.

    The system is defined by:
    - Three spins: S1, S2, S3, each taking values +1 or -1.
    - Energy: E = -J[S1*S2 + S1*S3 + S2*S3].
    - Partition function: Z = sum over all states of exp(-beta * E).

    The function returns "Correct" if the candidate answer matches the calculation,
    otherwise it returns a detailed explanation of the error.
    """
    
    # Define symbolic variables for building the mathematical expressions
    J, beta = sympy.symbols('J beta')

    # --- Step 1: Systematically derive the correct partition function ---

    # Generate all 2^3 = 8 possible microstates for the three spins
    spin_values = [-1, 1]
    all_states = list(product(spin_values, repeat=3))

    # Define the energy function as per the problem statement
    def calculate_energy(state, J_symbol):
        s1, s2, s3 = state
        return -J_symbol * (s1 * s2 + s1 * s3 + s2 * s3)

    # Calculate the energy for each of the 8 microstates
    energies = [calculate_energy(state, J) for state in all_states]

    # Group states by energy to find the degeneracy of each energy level.
    # A Counter is used to count occurrences of each unique energy value.
    # The result is a dictionary like {energy_level: degeneracy_count}.
    energy_degeneracies = Counter(energies)

    # Construct the partition function Z = sum(g_i * exp(-beta * E_i))
    # where g_i is the degeneracy of energy level E_i.
    calculated_Z = sympy.sympify(0)
    for energy, degeneracy in energy_degeneracies.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # --- Step 2: Define the candidate answer to be checked ---
    
    # The provided final answer is C: Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    candidate_answer_expr = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # --- Step 3: Compare the derived result with the candidate answer ---

    # sympy.simplify(A - B) == 0 is a robust method to check for symbolic equality.
    # It simplifies the difference, and if the result is zero, the expressions are equivalent.
    if sympy.simplify(calculated_Z - candidate_answer_expr) == 0:
        return "Correct"
    else:
        # If they don't match, construct a detailed explanation for the discrepancy.
        reason = "Incorrect. The candidate answer does not match the derived partition function.\n\n"
        reason += "REASONING:\n"
        reason += "1. The system has 3 spins, each taking values +1 or -1, resulting in 2^3 = 8 total microstates.\n"
        reason += "2. The energy for each state is calculated using E = -J[S1*S2 + S1*S3 + S2*S3].\n"
        reason += "3. By enumerating all 8 states, the code finds the following energy levels and their degeneracies:\n"
        
        # Sort items for a consistent and readable output
        sorted_degeneracies = sorted(energy_degeneracies.items(), key=lambda item: str(item[0]))
        for energy, degeneracy in sorted_degeneracies:
            reason += f"   - Energy Level E = {energy}: Found in {degeneracy} states.\n"
        
        reason += "\n4. The partition function is Z = sum(g_i * exp(-beta * E_i)). Based on the degeneracies above, the correct expression is constructed.\n"
        reason += f"\nDERIVED PARTITION FUNCTION:\n   Z = {calculated_Z}\n"
        reason += f"\nCANDIDATE ANSWER:\n   Z = {candidate_answer_expr}\n"
        reason += "\nCONCLUSION: The candidate answer's expression does not match the one derived from first principles."
        
        return reason

# To check the answer, you can run the following line:
# print(check_partition_function_correctness())