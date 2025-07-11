import math

def calculate_dimension():
    """
    Calculates the Minkowski–Bouligand dimension for the described fractal piano.
    The dimension D is given by the formula: D = log(45) / log(14),
    which can be expanded to D = (log(5) + log(9)) / log(14).
    """
    # The numbers in the derived formula
    num1 = 5  # Number of self-similar copies
    num2 = 9  # Ratio of the scaling factors (sx/sy)
    num3 = 14 # The scaling factor in the x-direction (sx = 1/rx)

    # Calculate the numerator and denominator of the dimension formula
    numerator = math.log(num1) + math.log(num2)
    denominator = math.log(num3)

    # Calculate the dimension
    dimension = numerator / denominator

    # Print the explanation and the result
    print("The Minkowski–Bouligand dimension (D) is calculated using the formula derived from box-counting for self-affine fractals.")
    print(f"The formula is: D = (log({num1}) + log({num2})) / log({num3})")
    print(f"D = log({num1 * num2}) / log({num3})")
    print(f"D = {numerator:.4f} / {denominator:.4f}")
    print(f"The calculated dimension is: {dimension:.4f}")

calculate_dimension()