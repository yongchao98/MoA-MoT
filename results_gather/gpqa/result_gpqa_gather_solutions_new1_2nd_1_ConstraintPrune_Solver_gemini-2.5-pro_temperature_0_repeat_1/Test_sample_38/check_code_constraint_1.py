import numpy as np
import sympy as sp
import re

def check_correctness():
    """
    This function verifies the solution to the quantum mechanics problem.
    It follows these steps:
    1. Defines the problem symbolically using sympy.
    2. Calculates the probability density |ψ(x)|² from the given wave function.
    3. Sets up the normalization integral: ∫|ψ(x)|² dx = 1 over the domain [1, 3].
    4. Solves the integral equation to find the exact value of 'a'.
    5. Calculates the numerical value of 'a'.
    6. Compares this value to the given options to find the closest match.
    7. Checks if the provided answer from the LLM corresponds to this correct match.
    """
    try:
        # Step 1: Define the problem symbolically
        x, a = sp.symbols('x a', real=True, positive=True)

        # The wave function is ψ(x) = (a / sqrt(1 + x)) - 0.5*i
        # Real part: Re(ψ) = a / sqrt(1 + x)
        # Imaginary part: Im(ψ) = -0.5

        # Step 2: Calculate the probability density |ψ(x)|² = Re(ψ)² + Im(ψ)²
        prob_density = (a / sp.sqrt(1 + x))**2 + (-0.5)**2
        # This simplifies to a**2 / (1 + x) + 0.25

        # Step 3: Set up and solve the normalization integral over the domain [1, 3]
        # The total probability must be 1.
        # ∫[1 to 3] (a**2 / (1 + x) + 0.25) dx = 1
        total_probability_integral = sp.integrate(prob_density, (x, 1, 3))

        # Step 4: Solve the equation for 'a'
        # The integral evaluates to a**2 * ln(2) + 0.5
        # We set this equal to 1 and solve for 'a'
        equation = sp.Eq(total_probability_integral, 1)
        solutions = sp.solve(equation, a)

        # Since 'a' is defined as positive, we expect one solution
        if not solutions:
            return "Checker Error: Failed to solve for 'a' symbolically."
        
        a_exact_value = solutions[0]
        a_numerical_value = a_exact_value.evalf()

        # Step 5: Check which option is closest to the calculated value
        # The options as listed in the question prompt
        options = {
            "A": 0.85,
            "B": 0.35,
            "C": 0.6,
            "D": 1.1
        }

        closest_option_label = None
        min_difference = float('inf')

        for label, value in options.items():
            difference = abs(value - a_numerical_value)
            if difference < min_difference:
                min_difference = difference
                closest_option_label = label
        
        correct_value = options[closest_option_label]

        # Step 6: Parse the provided LLM answer
        llm_final_answer = "<<<A>>>"
        match = re.search(r'<<<([A-D])>>>', llm_final_answer)
        if not match:
            return f"Checker Error: Could not parse the provided answer format: {llm_final_answer}"
            
        llm_answer_label = match.group(1)
        llm_answer_value = options.get(llm_answer_label)

        if llm_answer_value is None:
            return f"Checker Error: The provided answer label '{llm_answer_label}' is not a valid option."

        # Step 7: Compare the LLM's answer with the calculated correct answer
        if llm_answer_label == closest_option_label:
            return "Correct"
        else:
            return (f"Incorrect. The calculation shows that the value of 'a' must be approximately {a_numerical_value:.4f}. "
                    f"The closest option is {correct_value} (Option {closest_option_label}). "
                    f"The provided answer was {llm_answer_value} (Option {llm_answer_label}).")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)