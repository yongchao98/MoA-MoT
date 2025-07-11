import math

def predict_product_ratio():
    """
    This function predicts the ratio of isomers A and B based on the number of
    pi electrons in their respective formation pathways.
    """
    # The problem describes a thermal reaction of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene,
    # which produces a cis-isomer A and a trans-isomer B.
    # This is explained by two competing electrocyclic reactions:

    # Product A (cis-isomer) is formed via a 6-pi electron disrotatory cyclization.
    pi_electrons_A = 6

    # Product B (trans-isomer) is formed via an 8-pi electron conrotatory cyclization.
    pi_electrons_B = 8

    # To predict the ratio of A to B, we can use a simple model where the reaction rate
    # is inversely proportional to the number of pi electrons involved. This is a proxy
    # for the activation energy barrier, as the larger 8-pi system requires a more
    # sterically strained conformation to react.
    # So, Ratio(A/B) is proportional to pi_electrons_B / pi_electrons_A.

    # The numerator of our ratio corresponds to the number of pi electrons for product B.
    ratio_numerator = pi_electrons_B
    # The denominator of our ratio corresponds to the number of pi electrons for product A.
    ratio_denominator = pi_electrons_A

    # To express this as the simplest integer ratio, we can divide both numbers
    # by their greatest common divisor (GCD).
    common_divisor = math.gcd(ratio_numerator, ratio_denominator)

    final_ratio_A = ratio_numerator // common_divisor
    final_ratio_B = ratio_denominator // common_divisor

    print("Step 1: Identify the competing pathways and the number of pi electrons for each.")
    print(f"Pathway to A (cis-isomer): 6-pi electron cyclization. N_A = {pi_electrons_A}")
    print(f"Pathway to B (trans-isomer): 8-pi electron cyclization. N_B = {pi_electrons_B}")
    print("\nStep 2: Propose a model to estimate the product ratio.")
    print("Model: The ratio of products A:B is assumed to be N_B : N_A.")
    print("\nStep 3: Calculate the ratio.")
    print(f"The predicted ratio of A : B is {ratio_numerator} : {ratio_denominator}.")
    print(f"Simplifying this gives the final ratio: {final_ratio_A} : {final_ratio_B}.")

predict_product_ratio()
<<<4:3>>>