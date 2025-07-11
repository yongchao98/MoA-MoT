import math

def calculate_l(d, lambda_val):
    """
    Calculates the value of l(d, lambda) based on the derived analytical formula.

    The derivation is as follows:
    1. The probability density p_{d, lambda}(x) is found to be proportional to
       exp(-arccos(x^T . mu)^2 / (2 * lambda)), where mu = 1_d / sqrt(d).
    2. l(d, lambda) = ln(p(x1)/p(x2)) simplifies to
       (arccos(x2^T . mu)^2 - arccos(x1^T . mu)^2) / (2 * lambda).
    3. We compute the dot products:
       x1^T . mu = (e1+e2+e3)/sqrt(3) . 1_d/sqrt(d) = 3 / sqrt(3d) = sqrt(3/d).
       x2^T . mu = (e3+e4)/sqrt(2) . 1_d/sqrt(d) = 2 / sqrt(2d) = sqrt(2/d).
    4. This gives the final formula implemented below.
    """
    if d < 4:
        raise ValueError("d must be >= 4")
    if lambda_val < 1:
        raise ValueError("lambda must be >= 1")

    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lambda_val}")
    print("------------------------------------------")
    print("The final equation for l(d, lambda) is: (term2 - term1) / (2 * lambda)")
    print("")

    # Calculate the arguments for arccos
    cos_arg1 = math.sqrt(3 / d)
    cos_arg2 = math.sqrt(2 / d)

    # Calculate the squared arccos terms
    term1 = math.acos(cos_arg1)**2
    term2 = math.acos(cos_arg2)**2

    # Calculate the final result
    result = (term2 - term1) / (2 * lambda_val)

    # Output each number in the final equation
    print("Value of d: ", d)
    print("Value of lambda: ", lambda_val)
    print(f"cos_arg1 = sqrt(3/{d}) = {cos_arg1}")
    print(f"cos_arg2 = sqrt(2/{d}) = {cos_arg2}")
    print(f"term1 = arccos(sqrt(3/{d}))^2 = {term1}")
    print(f"term2 = arccos(sqrt(2/{d}))^2 = {term2}")
    print(f"Final result l({d}, {lambda_val}) = ({term2} - {term1}) / (2 * {lambda_val}) = {result}")
    
    return result

if __name__ == '__main__':
    # As no specific values for d and lambda were provided,
    # we will use d=5 and lambda=2 as an example, which satisfy the constraints
    # d >= 4 and lambda >= 1.
    d_val = 5
    lambda_val = 2
    final_answer = calculate_l(d_val, lambda_val)
    # The final answer is printed within the function,
    # but we will also print it here to be explicit.
    # print(f"\nFinal numerical answer: {final_answer}")