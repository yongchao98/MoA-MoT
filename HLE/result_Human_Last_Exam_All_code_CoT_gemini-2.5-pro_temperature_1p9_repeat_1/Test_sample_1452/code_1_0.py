import re

def main():
    """
    This script finds and prints the famous equation associated with the physicist Max Born.
    """
    
    scientist_name = "Max Born"
    
    # Max Born's famous equation is the Born rule, which provides the probability
    # interpretation of the wave function in quantum mechanics.
    # The equation is: rho = |psi|^2
    # The number in this equation is 2.
    
    # We will construct the equation string using the number.
    eq_part_1 = r"\( \rho = |\psi|^"
    number = 2
    eq_part_2 = r" \)"
    
    # Final equation in LaTeX format
    final_equation = f"{eq_part_1}{number}{eq_part_2}"
    
    print(f"Name: {scientist_name}")
    print(f"Equation: {final_equation}")

if __name__ == "__main__":
    main()