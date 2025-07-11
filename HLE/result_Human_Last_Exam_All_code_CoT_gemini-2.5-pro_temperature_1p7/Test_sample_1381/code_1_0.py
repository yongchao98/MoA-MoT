import math

def calculate_max_equilibria(N):
    """
    Calculates the maximum number of equilibria for the given Lotka-Volterra system with N species.

    Args:
        N (int): The number of species in the system. Must be a non-negative integer.

    Returns:
        int: The maximum number of possible equilibria.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: The number of species (N) must be a non-negative integer.")
        return None
        
    # The maximum number of equilibria is 2^N.
    # We use the power operator ** for exponentiation.
    base = 2
    exponent = N
    result = base ** exponent
    
    print(f"For a system with N = {N} species:")
    print(f"The analysis shows the maximum number of possible equilibria is given by the equation: base^N")
    # Outputting each number in the final equation as requested
    print(f"Plugging in the numbers: {base}**{exponent} = {result}")
    return result

# --- Example Usage ---
# We consider a system with N=10 species as an example.
# You can change this value to explore other cases.
number_of_species = 10
calculate_max_equilibria(number_of_species)
