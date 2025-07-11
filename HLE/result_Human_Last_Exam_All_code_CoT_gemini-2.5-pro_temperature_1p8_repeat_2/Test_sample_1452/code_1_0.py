import collections

def print_equations():
    """
    Prints a list of famous scientists and their equations, including Max Born's.
    """
    # Using an OrderedDict to maintain the order as provided in the prompt
    scientist_equations = collections.OrderedDict([
        ("Ernst Abbe", r"d = \frac{\lambda}{2 \, \text{n sin} \, \alpha}"),
        ("Ludwig Boltzmann", r"S = k \log W"),
        ("Erwin Schrödinger", r"i\hbar \frac{\partial}{\partial t} \Psi = \hat{H} \Psi"),
        ("Julian Schwinger", r"\alpha / 2 \pi"),
        ("Otto Hahn", r"^{92}\text{U} + ^1_0\text{n} \rightarrow"),
        ("P. A. M. Dirac", r"i \gamma^\mu \partial_\mu \psi = m \psi"),
        ("Stephen Hawking", r"T = \frac{\hbar c^3}{8 \pi G M k}"),
        ("Ferdinand von Lindemann", r"\pi"),
        ("Max Born", r"\int_{-\infty}^{\infty} |\psi(x,t)|^2 dx = 1")
    ])

    for name, equation in scientist_equations.items():
        print(f"Name: {name}")
        # The raw LaTeX string contains numbers like 2 and 1 which will be printed.
        print(f"Equation: \\( {equation} \\)")
        print()

print_equations()