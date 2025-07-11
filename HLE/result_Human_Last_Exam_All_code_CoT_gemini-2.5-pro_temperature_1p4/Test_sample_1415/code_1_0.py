import math

def solve_integral_problem():
    """
    This function solves the definite integral problem and reports the final value as requested.
    """
    # The integral to evaluate is ∫ (F₃'(x)/F₃(x)) dx from 0 to 1.
    # The antiderivative of f'(x)/f(x) is ln|f(x)|.
    # By the Fundamental Theorem of Calculus, the definite integral is ln|F₃(1)| - ln|F₃(0)|.
    # F₃(x) = ln(1 + e^(sin(x))). On [0, 1], F₃(x) > 0, so we can drop the absolute value.
    # The value of the integral is V = ln(F₃(1)) - ln(F₃(0)).

    # Step 1: Calculate F₃(0)
    # F₃(0) = ln(1 + e^(sin(0))) = ln(1 + e^0) = ln(1 + 1) = ln(2)
    f3_0 = math.log(2)

    # Step 2: Calculate F₃(1)
    # F₃(1) = ln(1 + e^(sin(1)))
    f3_1 = math.log(1 + math.exp(math.sin(1)))

    # Step 3: Calculate the value of the integral, V
    # V = ln(F₃(1)) - ln(F₃(0))
    V = math.log(f3_1) - math.log(f3_0)

    # Step 4: Calculate 10000 * V and find the closest integer
    final_result = 10000 * V
    closest_integer = round(final_result)
    
    # Print the breakdown of the calculation
    print("The integral is V = ln(F₃(1)) - ln(F₃(0)), where F₃(x) = ln(1 + e^(sin(x))).")
    print("\nFirst, we evaluate the terms:")
    print(f"F₃(0) = ln(2) = {f3_0}")
    print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_1}")
    
    print("\nPlugging these into the equation for V:")
    # Printing the numbers that form the final equation
    print(f"V = ln({f3_1}) - ln({f3_0})")
    print(f"The value of the integral is V = {V}")

    print("\nFinally, we calculate 10000 * V and find the closest integer:")
    print(f"10000 * V = {final_result}")
    print(f"The closest integer is: {closest_integer}")
    
    return closest_integer

if __name__ == '__main__':
    solve_integral_problem()
