import sympy
from collections import defaultdict

def check_partition_function_correctness():
    """
    This function programmatically calculates the partition function for the given 3-spin system
    and compares it to the expression provided in the final answer.
    """
    # Define symbolic variables for the calculation
    J, beta = sympy.symbols('J beta')

    # The system has three spins, each can be +1 or -1
    spin_values = [1, -1]
    
    # Use a dictionary to store the unique energy levels and count their degeneracies
    energy_degeneracy = defaultdict(int)

    # 1. Iterate through all 2^3 = 8 possible microstates of the system
    num_states = 0
    for s1 in spin_values:
        for s2 in spin_values:
            for s3 in spin_values:
                num_states += 1
                # 2. Calculate the energy for the current state using the given formula
                # E = -J[S1S2 + S1S3 + S2S3]
                energy = -J * (s1*s2 + s1*s3 + s2*s3)
                
                # 3. Increment the count (degeneracy) for this energy level
                energy_degeneracy[energy] += 1

    # Sanity check: ensure all 8 states were processed
    if num_states != 8:
        return f"Incorrect: The code failed to iterate through all 8 states. Found {num_states} states."

    # 4. Construct the partition function Z by summing the Boltzmann factors for each unique energy level
    # Z = Σ_levels g(E) * e^(-βE)
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # 5. Define the expression from the selected answer (B)
    # The final answer is Z = 2e^(3Jβ) + 6e^(-Jβ)
    answer_expression = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 6. Compare the programmatically calculated expression with the answer's expression.
    # sympy.simplify(expr1 - expr2) == 0 is a robust way to check for symbolic equality.
    if sympy.simplify(calculated_Z - answer_expression) == 0:
        return "Correct"
    else:
        # If they don't match, provide a detailed reason.
        reason = "Incorrect: The calculated partition function does not match the expression from the provided answer.\n"
        
        # Extract calculated degeneracies for clarity
        calculated_degeneracies = {str(k): v for k, v in energy_degeneracy.items()}
        reason += f"Calculated energy levels and their degeneracies: {calculated_degeneracies}\n"
        
        # Show the expected degeneracies based on the correct answer
        expected_degeneracies = {'-3*J': 2, 'J': 6}
        reason += f"Expected energy levels and their degeneracies: {expected_degeneracies}\n"
        
        reason += f"The code calculated the partition function to be: Z = {calculated_Z}\n"
        reason += f"The expression from the answer is: Z = {answer_expression}"
        return reason

# Execute the check and print the result
result = check_partition_function_correctness()
print(result)