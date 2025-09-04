import sympy
from itertools import product

def check_correctness_of_answer_A():
    """
    This function programmatically calculates the partition function for the given 
    three-spin system and verifies if the provided answer 'A' is correct.

    The system has three spins (S1, S2, S3), each taking values +1 or -1.
    The energy is given by E = -J[S1*S2 + S1*S3 + S2*S3].
    The partition function is Z = sum(exp(-beta * E)) over all possible states.
    
    To run this code, the 'sympy' library must be installed (`pip install sympy`).
    """
    
    # 1. Define symbolic variables for the constants J and beta.
    # Using a symbolic math library allows for exact comparison of expressions.
    J, beta = sympy.symbols('J beta')

    # 2. Define the possible values for a single spin.
    spin_values = [+1, -1]

    # 3. Generate all 2^3 = 8 possible microstates for the three-spin system.
    # itertools.product is perfect for generating these combinations.
    all_microstates = list(product(spin_values, repeat=3))
    
    # 4. Calculate the partition function by summing the Boltzmann factors for each state.
    Z_calculated = 0
    
    # This dictionary will store the energy levels and their degeneracies (counts)
    # for a more detailed analysis.
    energy_degeneracy = {}

    for state in all_microstates:
        s1, s2, s3 = state
        
        # Calculate the energy for the current microstate using the given formula.
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        
        # Add the Boltzmann factor exp(-beta*E) to the total sum for Z.
        Z_calculated += sympy.exp(-beta * energy)

        # Count the degeneracy of each energy level.
        energy_degeneracy[energy] = energy_degeneracy.get(energy, 0) + 1

    # 5. Simplify the resulting symbolic expression for the partition function.
    # The summation results in terms that can be grouped, e.g., exp(x) + exp(x) = 2*exp(x).
    Z_calculated_simplified = sympy.simplify(Z_calculated)

    # 6. Define the expression from the given answer (Option A).
    # Answer A is: Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    Z_from_answer_A = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 7. Compare the calculated result with the answer's expression.
    # If the expressions are mathematically equivalent, their difference will simplify to 0.
    if sympy.simplify(Z_calculated_simplified - Z_from_answer_A) == 0:
        # The calculated Z matches the answer's Z.
        return "Correct"
    else:
        # The calculated Z does not match the answer's Z.
        # Provide a detailed reason for the discrepancy.
        reason = (
            "The answer is incorrect.\n"
            f"Reasoning: The calculation was performed by summing the Boltzmann factor e^(-beta*E) over all 8 possible spin configurations.\n"
            f"The distinct energy levels and their degeneracies were found to be: {energy_degeneracy}.\n"
            f"This leads to a correctly calculated partition function: Z = {Z_calculated_simplified}.\n"
            f"The provided answer (A) corresponds to Z = {Z_from_answer_A}, which is not equivalent to the calculated result."
        )
        return reason

# You can run the function to get the result.
# print(check_correctness_of_answer_A())