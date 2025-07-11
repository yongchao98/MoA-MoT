import math

def calculate_finetuning_factor():
    """
    Calculates the multiplicative factor change in the fine-tuning measure
    based on the change in the energy scale of a particle accelerator.
    """
    
    # Define the initial and final energy scales.
    # Current scale (Lambda_1) is 8 TeV.
    lambda_1_tev = 8.0
    
    # Proposed new scale (Lambda_2) is 1.1 PeV.
    lambda_2_pev = 1.1
    
    # To calculate the ratio, we need to convert both scales to the same unit.
    # We will use TeV. 1 PeV = 1000 TeV.
    lambda_2_tev = lambda_2_pev * 1000
    
    # As reasoned in the plan, the multiplicative factor is the ratio of the two scales.
    # This is because the fine-tuning measure is directly proportional to the scale Lambda.
    factor = lambda_2_tev / lambda_1_tev
    
    # The problem asks to output each number in the final equation.
    print("The multiplicative factor is calculated by the ratio of the new energy scale to the current energy scale:")
    print(f"Factor = {lambda_2_tev} TeV / {lambda_1_tev} TeV")
    
    # Print the final result, rounded to two decimal places.
    print(f"\nThe fine-tuning according to the measure changes by a multiplicative factor of: {factor:.2f}")

calculate_finetuning_factor()