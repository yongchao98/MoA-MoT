import sympy
import itertools

def check_partition_function_correctness():
    """
    This function programmatically calculates the partition function for the given 3-spin system
    and checks if the result matches the logic of the provided answer.

    The provided answer's reasoning correctly identifies the energy levels and their degeneracies,
    which leads to option B. This code verifies that calculation independently.
    """
    
    # 1. Define symbolic variables for J and beta to handle the expressions algebraically.
    J, beta = sympy.symbols('J beta')

    # 2. List all possible configurations for the three spins (S1, S2, S3).
    #    Each spin can be +1 or -1, so there are 2^3 = 8 configurations.
    spins = [+1, -1]
    configurations = list(itertools.product(spins, repeat=3))

    # 3. For each configuration, calculate the energy and count the degeneracy for each energy level.
    energy_degeneracy = {}
    energy_formula = lambda s1, s2, s3: -J * (s1*s2 + s1*s3 + s2*s3)

    for s1, s2, s3 in configurations:
        energy = energy_formula(s1, s2, s3)
        # Add the energy to the dictionary and increment its count (degeneracy).
        energy_degeneracy[energy] = energy_degeneracy.get(energy, 0) + 1

    # Expected energy levels and degeneracies from the correct calculation:
    # E = -3J (for states +++ and ---), degeneracy = 2
    # E = +J (for all 6 other states), degeneracy = 6
    
    # 4. Construct the partition function Z = sum(g_i * exp(-beta * E_i))
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # 5. Define the expression for option B, which is the result of the provided answer's logic.
    option_B = 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta)

    # 6. Check if the calculated partition function matches option B.
    #    sympy.simplify() ensures that two mathematically equivalent expressions are recognized as such.
    if sympy.simplify(calculated_Z - option_B) == 0:
        # The independent calculation confirms that the logic in the provided answer is sound
        # and leads to the correct result (Option B).
        return "Correct"
    else:
        # This block would execute if the provided answer's logic was flawed.
        return (f"Incorrect. The reasoning in the provided answer is flawed. "
                f"The correct partition function, calculated from first principles, is Z = {calculated_Z}. "
                f"This does not match the result implied by the provided answer (Option B).")

# Run the check and print the result.
result = check_partition_function_correctness()
print(result)