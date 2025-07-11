import math

def get_cn(n):
    """
    Calculates the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram in the virial expansion.

    Args:
        n (int): The number of nodes (particles) in the diagram.

    Returns:
        float: The value of the prefactor c_n.
    """
    if n < 2:
        raise ValueError("n must be 2 or greater for virial coefficients.")
    
    if n == 2:
        # Special case for n=2
        return -0.5
    else:
        # General formula for n >= 3
        return ((-1)**n) / n

# Print the prefactors for n from 2 to 8 as examples
print("The prefactor c_n for various n are:")
for n_val in range(2, 9):
    cn_value = get_cn(n_val)
    print(f"c_{n_val} = {cn_value}")

# The user asked for a single formula. Let's summarize it as a final answer.
# For n=2, c_2 = -1/2
# For n>=3, c_n = (-1)^n / n
# Example values from the code:
# c_2 = -0.5
# c_3 = -0.333...
# c_4 = 0.25
# c_5 = -0.2
# c_6 = 0.166...
# c_7 = -0.142...
# c_8 = 0.125
# Since the format requested is <<<answer content>>>, I will provide the expression.
# The question "Determine the system-independent prefactor cn" implies finding the formula for it.
final_answer_expression = "c_n = -1/2 for n=2, and c_n = (-1)^n / n for n >= 3"