import sympy as sp

def display_charge_formulas():
    """
    This function displays the formulas for total volume and surface charges
    for the given concentric cylindrical electrode problem, as presented in one of the answer choices.
    """
    # Define symbolic variables
    V, epsilon, pi, L, a, b = sp.symbols('V varepsilon pi L a b')
    
    # Formulas from Answer Choice B
    q_v = -4 * V * epsilon * pi * L / (1 - a**2/b**2)
    q_s_inner = 2 * pi * L * V * epsilon / (1 - a**2/b**2)
    q_s_outer = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2/b**2))
    
    # Print the formulas
    print("Formulas from the selected answer choice:")
    print(f"Total volume charge (q_v) = {sp.pretty(q_v)}")
    print(f"Total surface charge on inner electrode (q_s(r=a)) = {sp.pretty(q_s_inner)}")
    print(f"Total surface charge on outer electrode (q_s(r=b)) = {sp.pretty(q_s_outer)}")

if __name__ == "__main__":
    display_charge_formulas()