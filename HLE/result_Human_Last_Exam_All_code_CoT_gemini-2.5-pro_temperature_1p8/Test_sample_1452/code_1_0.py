import sympy

# The user wants to find the famous equation associated with Max Born.
# Max Born is famous for the Born rule, which gives the probability density
# of finding a particle at a given point in terms of its wavefunction Psi.
# The probability density is |Psi|^2.
# A common way to express this rule is through the normalization condition,
# which states that the total probability of finding the particle somewhere in all of space is 1.

# Define the symbols for the equation using LaTeX format.
# Although sympy can render LaTeX, for this task, we will just format the string.
scientist_name = "Max Born"
# LaTeX string for the Born rule (normalization condition)
equation_latex = r"\int_{-\infty}^{\infty} |\Psi(x)|^2 \, dx = 1"

# Print the result in the desired format.
print(f'Name: {scientist_name}')
print(f'Equation: {equation_latex}')