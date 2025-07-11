import math

def solve_integral():
    """
    This function solves the definite integral problem and reports the final value.
    """
    # Let F₃(x) = ln(1 + e^(sin(x)))
    # The integral is ∫ (F₃'(x)/F₃(x)) dx from 0 to 1.
    # The antiderivative of F₃'(x)/F₃(x) is ln(F₃(x)).
    # By the Fundamental Theorem of Calculus, the definite integral is ln(F₃(1)) - ln(F₃(0)).

    print("Step 1: Define the function F₃(x) and evaluate it at the limits of integration (x=1 and x=0).")
    
    # F₃(x) = ln(1 + e^(sin(x)))
    def F3(x):
        return math.log(1 + math.exp(math.sin(x)))

    # Evaluate at x = 1
    F3_at_1 = F3(1)
    
    # Evaluate at x = 0
    F3_at_0 = F3(0)
    
    print(f"The upper limit value is F₃(1) = ln(1 + e^(sin(1))) = {F3_at_1}")
    print(f"The lower limit value is F₃(0) = ln(1 + e^(sin(0))) = {F3_at_0}\n")
    
    # Since F₃(x) > 0 for x in [0, 1], we can drop the absolute value in ln|F₃(x)|.
    
    print("Step 2: Calculate the value of the integral V = ln(F₃(1)) - ln(F₃(0)).")
    
    # V = ln(F₃(1)) - ln(F₃(0))
    V = math.log(F3_at_1) - math.log(F3_at_0)
    
    # Print the equation with the calculated numbers
    print(f"V = ln({F3_at_1}) - ln({F3_at_0})")
    print(f"V = {V}\n")

    print("Step 3: Calculate the final result as the closest integer to 10000 * V.")
    
    final_value = 10000 * V
    closest_integer = round(final_value)
    
    print(f"10000 * V = 10000 * {V} = {final_value}")
    print(f"The closest integer is: {closest_integer}")

solve_integral()
<<<5487>>>