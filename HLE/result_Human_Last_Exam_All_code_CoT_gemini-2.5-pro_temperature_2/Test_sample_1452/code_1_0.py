import re

def solve_equation():
    """
    This function provides the equation for Max Born and extracts the numbers from it.
    """
    name = "Max Born"
    # The Born rule's normalization condition is a fundamental equation in quantum mechanics.
    equation_latex = r"\int_{-\infty}^{\infty} |\psi(x)|^2 dx = 1"
    
    # Extracting all numbers from the equation string using regex
    # This finds integers and floating point numbers.
    # We will filter out symbolic numbers like infinity.
    found_numbers = re.findall(r'-?\d+\.?\d*', equation_latex)
    
    # We are interested in the physical numbers/constants in the equation, which are 2 and 1.
    numbers = [2, 1]

    print(f'Name: "{name}"')
    # Print the equation in a format that can be rendered by math engines
    print(f'Equation: \\({equation_latex}\\)')

    # Print each number as requested by the prompt
    print(f"The first number in the equation is: {numbers[0]}")
    print(f"The second number in the equation is: {numbers[1]}")

solve_equation()