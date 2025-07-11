import math

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor by which the fine-tuning measure changes.
    """
    # Define the initial and final energy scales
    Lambda1_TeV = 8.0
    Lambda2_PeV = 1.1

    # Convert both scales to a common unit, GeV
    # 1 TeV = 1e3 GeV
    # 1 PeV = 1e6 GeV
    Lambda1_GeV = Lambda1_TeV * 1e3
    Lambda2_GeV = Lambda2_PeV * 1e6

    # The core of the problem is realizing that the fine-tuning measure Delta
    # is directly proportional to the bare Higgs mass, which in turn is
    # directly proportional to the energy scale Lambda due to the quadratic
    # divergence of the top loop correction.
    # Therefore, the multiplicative factor is simply the ratio of the two scales.
    
    print("The multiplicative factor is the ratio of the new energy scale to the old one.")
    print("The equation is: Factor = Lambda_2 / Lambda_1")
    print(f"Plugging in the values in GeV: Factor = {Lambda2_GeV:.0f} / {Lambda1_GeV:.0f}")

    # Calculate the ratio
    factor = Lambda2_GeV / Lambda1_GeV

    # Print the final result, rounded to two decimal places
    print(f"\nThe multiplicative factor is: {factor:.2f}")

calculate_fine_tuning_factor()
<<<137.50>>>