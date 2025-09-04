import sympy as sp
import numpy as np

def check_correctness():
    """
    This function verifies the solution to the quantum mechanics problem by:
    1. Symbolically deriving the expression for the normalization constant 'a'.
    2. Calculating the numerical value of 'a'.
    3. Comparing the calculated value to the given options to find the best fit.
    4. Checking if the best fit matches the provided answer.
    """
    
    # --- Model the problem constraints ---
    
    # The question provides the options:
    # A) 0.85, B) 0.35, C) 0.6, D) 1.1
    options = {
        "A": 0.85,
        "B": 0.35,
        "C": 0.6,
        "D": 1.1
    }
    
    # The final answer given by the LLM is <<<A>>>.
    llm_answer_label = "A"
    
    # --- Step 1: Symbolically derive the value of 'a' ---
    
    x, a = sp.symbols('x a', real=True, positive=True)
    
    # The wave function is ψ(x) = (a / sqrt(1 + x)) - 0.5*i
    # The probability density is |ψ(x)|² = (Real part)² + (Imaginary part)²
    prob_density = (a / sp.sqrt(1 + x))**2 + (sp.S(-1)/2)**2
    
    # The normalization condition states that the integral of the probability
    # density over the allowed region [1, 3] must be 1.
    integral_result = sp.integrate(prob_density, (x, 1, 3))
    
    # Set the integral equal to 1 to create the normalization equation.
    normalization_equation = sp.Eq(integral_result, 1)
    
    # Solve the equation for 'a'. We take the positive solution.
    solutions = sp.solve(normalization_equation, a)
    a_calculated_symbolic = [sol for sol in solutions if sol.is_positive][0]
    a_calculated_numerical = a_calculated_symbolic.evalf()

    # --- Step 2: Verify the LLM's answer ---
    
    # Find which of the given options is numerically closest to our calculated 'a'.
    closest_option_label = min(options.keys(), key=lambda k: abs(options[k] - a_calculated_numerical))
    
    # Check if the LLM's chosen answer label matches the label of the closest option.
    if llm_answer_label == closest_option_label:
        # As a final sanity check, let's confirm that the probability for the chosen 'a' is indeed close to 1.
        chosen_a_value = options[llm_answer_label]
        total_probability = integral_result.subs(a, chosen_a_value).evalf()
        
        if np.isclose(total_probability, 1.0, atol=1e-3):
            return "Correct"
        else:
            # This case would indicate an issue with the problem's options, but the logic is sound.
            return f"The chosen answer {llm_answer_label} is the closest, but it results in a total probability of {total_probability:.4f}, which is not sufficiently close to 1."
    else:
        return (f"Incorrect. The final answer was <<<{llm_answer_label}>>>, corresponding to a = {options[llm_answer_label]}. "
                f"However, the correct value for 'a' is calculated to be approximately {a_calculated_numerical:.4f}. "
                f"The closest option is {options[closest_option_label]} (Option {closest_option_label}). "
                f"Therefore, the correct answer should have been <<<{closest_option_label}>>>.")

# Run the check and print the result.
result = check_correctness()
print(result)