import sympy
import itertools

def check_answer():
    """
    This function calculates the partition function for the given spin system from first principles
    and compares it to the provided answer.
    """
    # Define symbolic variables for the calculation
    J, beta = sympy.symbols('J beta')

    # The final answer provided by the LLM is 'A', which corresponds to the expression:
    # Z = 2 * e^(3*J*beta) + 6 * e^(-J*beta)
    # Let's create a symbolic representation of this answer.
    expected_answer_expr = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # --- Start calculation from first principles ---

    # Each spin can be +1 or -1
    spins = [+1, -1]
    
    # Initialize the partition function Z as a symbolic expression
    calculated_Z = 0
    
    # Dictionary to store energy levels and their degeneracies for analysis
    energy_degeneracy = {}

    # There are 3 spins, so there are 2^3 = 8 possible states.
    # We iterate through all possible states (s1, s2, s3).
    all_states = list(itertools.product(spins, repeat=3))
    
    for s1, s2, s3 in all_states:
        # Calculate the energy for the current state using the given Hamiltonian:
        # E = -J[ S1*S2 + S1*S3 + S2*S3 ]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        
        # The partition function is the sum of Boltzmann factors over all states:
        # Z = sum(exp(-beta * E))
        calculated_Z += sympy.exp(-beta * energy)
        
        # Keep track of the energy levels and their degeneracies
        if energy in energy_degeneracy:
            energy_degeneracy[energy] += 1
        else:
            energy_degeneracy[energy] = 1
            
    # Simplify the calculated expression to group terms with the same exponent.
    # sympy.simplify() will combine the terms into a canonical form.
    simplified_calculated_Z = sympy.simplify(calculated_Z)

    # --- Verification Step ---

    # Compare the calculated partition function with the expected answer.
    # sympy's equality operator (==) checks for symbolic equivalence.
    if simplified_calculated_Z == expected_answer_expr:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason = (
            f"The provided answer is incorrect.\n"
            f"The calculated partition function is: Z = {simplified_calculated_Z}\n"
            f"The provided answer expression is: Z = {expected_answer_expr}\n\n"
            "Analysis of the calculation from first principles reveals the following energy levels and degeneracies:\n"
        )
        
        # Sort the energies for a clear output
        sorted_energies = sorted(energy_degeneracy.items(), key=lambda item: item[0].as_coeff_Mul()[0])

        for energy, degeneracy in sorted_energies:
            reason += f"- Energy Level E = {energy}: Found {degeneracy} states (degeneracy = {degeneracy}).\n"
            
        reason += "\nThis leads to the correct partition function Z = 2*exp(3*J*beta) + 6*exp(-J*beta), which does not match the provided answer's expression."
        
        return reason

# Execute the check and print the result
print(check_answer())