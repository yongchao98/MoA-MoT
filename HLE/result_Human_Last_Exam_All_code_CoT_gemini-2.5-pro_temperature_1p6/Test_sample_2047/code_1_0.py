import math

def calculate_l(d, lambda_val):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    The formula for l(d, lambda) is:
    l(d, lambda) = (1/(2*lambda)) * [ arccos(x2.mu)^2 - arccos(x1.mu)^2 ]
    where:
    x1.mu = sqrt(3)/sqrt(d)
    x2.mu = sqrt(2)/sqrt(d)
    """

    if d < 4:
        print("Error: d must be an integer greater than or equal to 4.")
        return
    if lambda_val < 1:
        print("Error: lambda must be a float greater than or equal to 1.")
        return

    # Calculate the components of the formula
    dot_x1_mu = math.sqrt(3) / math.sqrt(d)
    dot_x2_mu = math.sqrt(2) / math.sqrt(d)

    # Ensure arguments for arccos are in the valid range [-1, 1]
    # This is guaranteed by the constraints on d.
    
    acos_x1_mu_sq = math.acos(dot_x1_mu) ** 2
    acos_x2_mu_sq = math.acos(dot_x2_mu) ** 2
    
    two_lambda = 2 * lambda_val

    # Calculate the final result
    result = (1 / two_lambda) * (acos_x2_mu_sq - acos_x1_mu_sq)

    # Print out each number in the final equation as requested
    print(f"For d = {d} and lambda = {lambda_val}:")
    print(f"l(d, lambda) = (1 / (2 * {lambda_val})) * [ (arccos(sqrt(2)/sqrt({d})))^2 - (arccos(sqrt(3)/sqrt({d})))^2 ]")
    print(f"             = (1 / {two_lambda}) * [ (arccos({dot_x2_mu}))^2 - (arccos({dot_x1_mu}))^2 ]")
    print(f"             = (1 / {two_lambda}) * [ {acos_x2_mu_sq} - {acos_x1_mu_sq} ]")
    print(f"Result: {result}")


if __name__ == '__main__':
    # Use example values d=4 and lambda=1.0 as they are not specified in the prompt.
    d_val = 4
    lambda_val = 1.0
    calculate_l(d_val, lambda_val)
