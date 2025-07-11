import math

def solve_integral():
    """
    Solves the definite integral problem step-by-step.
    The integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 simplifies to ln(F₃(1)) - ln(F₃(0)).
    This function calculates this value, multiplies by 10000, and rounds to the nearest integer.
    """
    
    # Define the chain of functions
    def F0(x):
        return x

    def F1(x):
        return math.sin(F0(x))

    def F2(x):
        return math.exp(F1(x))

    def F3(x):
        return math.log(1 + F2(x))

    # Evaluate F₃(x) at the upper and lower bounds of the integral
    F3_at_1 = F3(1)
    F3_at_0 = F3(0)

    # The value of the integral V = ln(F₃(1)) - ln(F₃(0))
    V = math.log(F3_at_1) - math.log(F3_at_0)

    # The final result is the closest integer to 10000 * V
    final_answer = round(10000 * V)

    print("The integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 evaluates to V = ln(F₃(1)) - ln(F₃(0)).")
    print("\n--- Calculation Breakdown ---")
    
    # Print the equation parts
    print(f"F₃(1) = ln(1 + e^(sin(1))) ≈ {F3_at_1}")
    print(f"F₃(0) = ln(2) ≈ {F3_at_0}")
    print(f"\nTherefore, V = ln({F3_at_1}) - ln({F3_at_0})")
    print(f"V ≈ {V}")

    print("\n--- Final Answer Calculation ---")
    print(f"The required value is the closest integer to 10000 * V.")
    print(f"10000 * V ≈ {10000 * V}")
    print(f"The closest integer is: {final_answer}")


solve_integral()
<<<5487>>>