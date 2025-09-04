import sympy
import itertools

def check_answer():
    """
    This function programmatically calculates the partition function for a system of 
    three spins and checks if it matches the provided answer.
    """
    # Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # Each spin can be +1 or -1
    spin_values = [+1, -1]

    # Generate all 2^3 = 8 possible configurations for the three spins
    all_configurations = list(itertools.product(spin_values, repeat=3))

    # Initialize the partition function Z to zero
    calculated_Z = 0

    # A dictionary to store energy levels and their degeneracies (counts)
    energy_degeneracy = {}

    # Loop over all possible configurations
    for config in all_configurations:
        S1, S2, S3 = config
        
        # Calculate the energy for the current configuration using the given formula
        energy = -J * (S1*S2 + S1*S3 + S2*S3)
        
        # Add the Boltzmann factor e^(-beta*E) to the total partition function
        calculated_Z += sympy.exp(-beta * energy)
        
        # Count the degeneracy of each energy level
        energy_degeneracy[energy] = energy_degeneracy.get(energy, 0) + 1

    # The expression for the partition function from the provided answer (Option B)
    answer_Z = 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta)

    # To verify correctness, we simplify the difference between our calculated Z
    # and the answer's Z. If the result is 0, they are symbolically equivalent.
    if sympy.simplify(calculated_Z - answer_Z) == 0:
        return "Correct"
    else:
        # If they are not equivalent, provide a detailed reason for the mismatch.
        calculated_Z_simplified = sympy.collect(calculated_Z, sympy.exp(J*beta))
        
        reason = (
            "The answer is incorrect.\n"
            f"The provided answer is: Z = {answer_Z}\n"
            f"My calculation from first principles gives: Z = {calculated_Z_simplified}\n"
            "Reasoning:\n"
            f"There are {len(all_configurations)} total spin configurations.\n"
            "The calculated energy levels and their degeneracies are:\n"
        )
        for energy, degeneracy in energy_degeneracy.items():
            reason += f"  - Energy E = {energy}, Degeneracy g = {degeneracy}\n"
        
        reason += "The two expressions for Z are not equivalent."
        return reason

# Run the check and print the result.
result = check_answer()
print(result)