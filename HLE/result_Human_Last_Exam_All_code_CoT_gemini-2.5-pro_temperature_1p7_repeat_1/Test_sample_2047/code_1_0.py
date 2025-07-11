import math

def solve_l_d_lambda():
    """
    This function calculates the value of l(d, lambda) based on the derived formula.
    The values for d and lambda are set to the minimums from the function's domain,
    d >= 4 and lambda >= 1, for a specific numerical answer.
    """
    # Parameters from the function definition l(d, lambda)
    # d is an integer >= 4
    # lambda is a real number >= 1
    # Using example values d=4, lambda=1.0
    d = 4
    lmbda = 1.0

    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lmbda}")
    print("-" * 30)

    # Calculate intermediate terms
    # Term for x_2: arccos(sqrt(2/d))^2
    arg1_val = math.sqrt(2.0 / d)
    acos1_val = math.acos(arg1_val)
    term1_sq = acos1_val**2

    # Term for x_1: arccos(sqrt(3/d))^2
    arg2_val = math.sqrt(3.0 / d)
    acos2_val = math.acos(arg2_val)
    term2_sq = acos2_val**2

    # Denominator
    denominator = 2.0 * lmbda

    # Final result
    result = (term1_sq - term2_sq) / denominator
    
    # Outputting the numbers in the final equation as requested
    print("Final Equation:")
    print(f"l({d}, {lmbda}) = ( (arccos(sqrt(2/{d})))^2 - (arccos(sqrt(3/{d})))^2 ) / (2 * {lmbda})")
    print("\nStep-by-step calculation:")
    print(f"sqrt(2/d) = sqrt(2/{d}) = {arg1_val}")
    print(f"arccos(sqrt(2/d)) = {acos1_val}")
    print(f"(arccos(sqrt(2/d)))^2 = {term1_sq}\n")

    print(f"sqrt(3/d) = sqrt(3/{d}) = {arg2_val}")
    print(f"arccos(sqrt(3/d)) = {acos2_val}")
    print(f"(arccos(sqrt(3/d)))^2 = {term2_sq}\n")

    print(f"Numerator = {term1_sq} - {term2_sq} = {term1_sq - term2_sq}")
    print(f"Denominator = 2 * {lmbda} = {denominator}")
    
    print("\nFinal Result:")
    print(f"l({d}, {lmbda}) = {result}")

if __name__ == '__main__':
    solve_l_d_lambda()