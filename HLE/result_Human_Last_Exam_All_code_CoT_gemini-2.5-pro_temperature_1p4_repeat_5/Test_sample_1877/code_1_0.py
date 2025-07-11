import math

def solve_problem():
    """
    Calculates T(n,p) for n=3, p=7 and then finds ceil(T/12).
    """
    # Era B precision (p = 7 significant digits)
    p = 7

    # Divisor from the problem statement
    divisor = 12

    # A characteristic Lyapunov exponent for the chaotic regime of the logistic map is lambda = ln(2)
    lambda_val = math.log(2)
    
    # The time T is when an initial error of 10**(-p) grows to order 1
    # T = (p * ln(10)) / lambda
    ln10 = math.log(10)
    T = (p * ln10) / lambda_val

    # The final answer is ceil(T / 12)
    result = math.ceil(T / divisor)

    # Print the explanation and the final equation
    print("The formula to find the number of iterations T is: T = (p * ln(10)) / λ")
    print(f"Here, p = {p} (for 7-digit precision) and we use the characteristic Lyapunov exponent λ = ln(2).")
    print("\nCalculation steps:")
    print(f"p = {p}")
    print(f"ln(10) ≈ {ln10:.4f}")
    print(f"λ = ln(2) ≈ {lambda_val:.4f}")
    print(f"T ≈ ({p} * {ln10:.4f}) / {lambda_val:.4f} ≈ {T:.4f}")
    
    print("\nThe problem asks for ceil(T / 12). The final equation is:")
    # We output each number in the final equation as requested.
    print(f"ceil({T:.4f} / {divisor}) = {result}")


solve_problem()
