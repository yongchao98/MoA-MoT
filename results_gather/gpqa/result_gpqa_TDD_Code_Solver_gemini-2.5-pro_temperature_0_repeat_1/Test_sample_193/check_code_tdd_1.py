import sympy
import itertools

def check_partition_function_correctness():
    """
    This function programmatically calculates the partition function for the given 3-spin system
    and verifies if the expression from option C is correct.

    The system has three spins (S1, S2, S3), each taking values +1 or -1.
    The energy is given by E = -J[S1*S2 + S1*S3 + S2*S3].
    The partition function is Z = sum(exp(-beta * E)) over all possible states.
    """
    
    # Define symbolic variables for the constants J and beta
    J, beta = sympy.symbols('J beta')

    # Each spin can be in one of two states: +1 or -1
    spin_values = [1, -1]

    # Generate all 2^3 = 8 possible configurations of the three spins
    all_states = list(itertools.product(spin_values, repeat=3))

    # Initialize the partition function Z to zero
    calculated_Z = 0
    
    # For providing a detailed explanation, we can also track the energy levels and their degeneracies
    energy_degeneracy_map = {}

    # Iterate over every possible state of the system
    for state in all_states:
        s1, s2, s3 = state
        
        # Calculate the energy E for the current state using the given formula
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        
        # Add the Boltzmann factor, exp(-beta * E), for this state to the partition function
        calculated_Z += sympy.exp(-beta * energy)
        
        # Tally the degeneracy for the calculated energy level
        # We use simplify to ensure that expressions like -J*(-1) become J
        simplified_energy = sympy.simplify(energy)
        energy_degeneracy_map[simplified_energy] = energy_degeneracy_map.get(simplified_energy, 0) + 1

    # The expression from option C is Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    # Let's create a symbolic representation of this answer
    answer_C_expression = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # The most reliable way to check for symbolic equality is to simplify the difference
    # between the calculated expression and the target expression. If the result is 0, they are equivalent.
    if sympy.simplify(calculated_Z - answer_C_expression) == 0:
        return "Correct"
    else:
        # If the check fails, construct a detailed error message
        
        # Reconstruct Z from the calculated degeneracies to show the derivation
        derived_Z_from_degeneracy = sum(
            degeneracy * sympy.exp(-beta * energy) 
            for energy, degeneracy in energy_degeneracy_map.items()
        )

        reason = (
            f"Incorrect. The expression from option C does not match the derived partition function.\n"
            f"Constraint check failed for the expression: Z = {answer_C_expression}\n\n"
            f"DERIVATION DETAILS:\n"
            f"The system has the following energy levels and degeneracies (number of states with that energy):\n"
        )
        for energy, degeneracy in sorted(energy_degeneracy_map.items(), key=lambda item: str(item[0])):
             reason += f"  - Energy Level E = {energy}: Degeneracy g = {degeneracy}\n"
        
        reason += (
            f"\nBased on this, the correct partition function Z = sum(g * exp(-beta*E)) is:\n"
            f"Z = {sympy.simplify(derived_Z_from_degeneracy)}"
        )
        return reason

# Execute the checking function and print its output.
# This will return "Correct" if option C matches the calculation.
result = check_partition_function_correctness()
print(result)