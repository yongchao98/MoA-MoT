def print_equations():
    """
    Prints a list of famous scientists and their equations, including Max Born's.
    """
    scientists = [
        {"name": "Ernst Abbe", "equation": r"\[d = \frac{\lambda}{2 \, \text{n sin} \, \alpha}\]"},
        {"name": "Ludwig Boltzmann", "equation": r"\( S = k \log W \)"},
        {"name": "Erwin Schrödinger", "equation": r"\( i\hbar \frac{\partial}{\partial t} \Psi = \hat{H} \Psi \)"},
        {"name": "Julian Schwinger", "equation": r"\(\alpha / 2 \pi\)"},
        {"name": "Otto Hahn", "equation": r"\(^{92}\text{U} + ^1_0\text{n} \rightarrow\)"},
        {"name": "P. A. M. Dirac", "equation": r"\(i \gamma^\mu \partial_\mu \psi = m \psi\)"},
        {"name": "Stephen Hawking", "equation": r"\( T = \frac{\hbar c^3}{8 \pi G M k} \)"},
        {"name": "Ferdinand von Lindemann", "symbol": r"\(\pi\)"},
        {"name": "Max Born", "equation": r"\(|\Psi|^2\)"}
    ]

    for item in scientists:
        print(f"Name: {item['name']}")
        if "equation" in item:
            print(f"Equation: {item['equation']}")
        elif "symbol" in item:
            print(f"Symbol: {item['symbol']}")
        # Ensure a blank line between entries for readability, except for the last one.
        if item != scientists[-1]:
            print()

if __name__ == "__main__":
    print_equations()