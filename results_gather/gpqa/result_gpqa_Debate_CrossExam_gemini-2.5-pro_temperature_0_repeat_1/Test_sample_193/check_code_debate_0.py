import sympy
import itertools

def check_partition_function_correctness():
    """
    This function programmatically calculates the partition function for the given 3-spin system
    and compares it with the provided answer (Option D).
    """
    # Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # The system has 3 spins, each can be +1 or -1.
    # We generate all 2^3 = 8 possible microstates.
    spin_values = [1, -1]
    all_states = list(itertools.product(spin_values, repeat=3))

    # The energy of the system is given by E = -J[S1*S2 + S1*S3 + S2*S3].
    # The partition function is Z = sum over all states of exp(-beta * E).

    # We will find the unique energy levels and their degeneracies (counts).
    energy_degeneracy = {}
    for state in all_states:
        s1, s2, s3 = state
        
        # Calculate the energy for the current state
        energy = -J * (s1 * s2 + s1 * s3 + s2 * s3)
        
        # Store the energy and count its occurrences (degeneracy)
        energy_degeneracy[energy] = energy_degeneracy.get(energy, 0) + 1

    # Construct the calculated partition function from the energy levels and their degeneracies.
    Z_calculated = 0
    for energy, degeneracy in energy_degeneracy.items():
        Z_calculated += degeneracy * sympy.exp(-beta * energy)

    # The proposed answer from option D is Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    Z_proposed = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # To check for correctness, we verify if the calculated expression is equivalent
    # to the proposed one. The simplest way is to check if their difference simplifies to zero.
    if sympy.simplify(Z_calculated - Z_proposed) == 0:
        return "Correct"
    else:
        # If they are not equal, the answer is incorrect.
        # We provide the calculated degeneracies as the reason for the failure.
        expected_energy_degeneracy = {
            -3*J: 2,  # From the term 2*exp(-beta*(-3J))
            J: 6      # From the term 6*exp(-beta*(J))
        }
        
        # Format the dictionaries for a clear error message
        calculated_str = {str(k): v for k, v in energy_degeneracy.items()}
        expected_str = {str(k): v for k, v in expected_energy_degeneracy.items()}
        
        reason = (
            "The answer is incorrect because the degeneracies of the energy levels are wrong.\n"
            f"Calculated (Energy: Degeneracy): {calculated_str}\n"
            f"Implied by Answer (Energy: Degeneracy): {expected_str}"
        )
        return reason

# Execute the check and print the result
result = check_partition_function_correctness()
print(result)