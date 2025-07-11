import sys

def display_greens_function():
    """
    Displays the functional form of the bare Green's function G_0
    and its dependence on the single-particle energy eigenvalue epsilon_k.
    """

    print("In the Feynman path integral formalism, the time-ordered bare (non-interacting) Green's function, G_0(k, omega), has the following functional form:")
    print("")
    print("         1")
    print("G_0 = ----------------------------------------")
    print("      (omega - (epsilon_k - mu)) + i*eta*sgn(epsilon_k - mu)")
    print("")
    print("This shows that G_0 is inversely proportional to the frequency (omega) shifted by the single-particle energy (epsilon_k).")
    print("\n--- Components of the Equation ---")

    # The prompt asks to "output each number in the final equation!".
    # Since the general form has no numbers, we will output each symbolic component.
    print(f"G_0: The bare Green's function, which is a function of momentum k and frequency omega.")
    print(f"omega: The frequency of the particle/excitation.")
    print(f"epsilon_k: The single-particle energy eigenvalue for momentum k.")
    print(f"mu: The chemical potential, which sets the zero point for energy.")
    print(f"i: The imaginary unit.")
    print(f"eta: A small positive infinitesimal (eta -> 0+) that ensures causality by shifting the poles of the function off the real axis.")
    print(f"sgn(x): The signum function, which is +1 for x > 0 and -1 for x < 0. It determines how the poles are shifted for particles vs. holes.")

if __name__ == '__main__':
    display_greens_function()