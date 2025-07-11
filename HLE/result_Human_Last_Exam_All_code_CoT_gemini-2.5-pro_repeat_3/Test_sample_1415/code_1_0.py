import math

def solve_integral():
    """
    Solves the definite integral and finds the closest integer to 10000*V.
    """
    # The integral to evaluate is ∫ (F₃'(x)/F₃(x)) dx from 0 to 1.
    # The antiderivative of f'(x)/f(x) is ln(f(x)).
    # So the definite integral is ln(F₃(1)) - ln(F₃(0)).

    # Step 1: Calculate F₃(0)
    # F₀(0) = 0
    # F₁(0) = sin(0) = 0
    # F₂(0) = e⁰ = 1
    # F₃(0) = ln(1 + 1) = ln(2)
    val_at_0 = math.log(2)

    # Step 2: Calculate F₃(1)
    # F₀(1) = 1
    # F₁(1) = sin(1)
    # F₂(1) = e^(sin(1))
    # F₃(1) = ln(1 + e^(sin(1)))
    val_at_1 = math.log(1 + math.exp(math.sin(1)))

    # Step 3: Calculate the value of the integral V
    # V = ln(F₃(1)) - ln(F₃(0))
    V = math.log(val_at_1) - math.log(val_at_0)

    # Step 4: Calculate 10000 * V and find the closest integer
    final_value = 10000 * V
    closest_integer = round(final_value)

    # Print the breakdown of the final equation as requested
    print("The definite integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 evaluates to ln(F₃(1)) - ln(F₃(0)).")
    print(f"The value of F₃(1) is: {val_at_1}")
    print(f"The value of F₃(0) is: {val_at_0}")
    print(f"The final value of the integral V = ln({val_at_1}) - ln({val_at_0}) is: {V}")
    print(f"The value of 10000 * V is: {final_value}")
    print(f"The closest integer to 10000 * V is: {closest_integer}")
    
    return closest_integer

# Run the solver and store the final answer
final_answer = solve_integral()
# The final answer is wrapped as requested
# print(f"\n<<<{final_answer}>>>")