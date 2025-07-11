import math

def calculate_l(n: int) -> float:
    """
    Calculates the value of l(n) for a given integer n >= 5.

    The function implements the derived analytical formula:
    l(n) = (2 * (n**2 + 1 - (2*n - 1) * sqrt(n**2 - n + 1))) / n**2
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    term_sqrt = math.sqrt(n**2 - n + 1)
    numerator = 2 * (n**2 + 1 - (2*n - 1) * term_sqrt)
    denominator = n**2
    
    return numerator / denominator

def print_equation_components():
    """
    Prints the components of the final equation for l(n) as requested.
    The equation is of the form:
    (C1*n**2 + C2 - (C3*n + C4)*sqrt(C5*n**2 + C6*n + C7)) / (C8*n**2)
    """
    # The derived formula is (2*n**2 + 2 - 2*(2*n - 1)*sqrt(n**2 - n + 1)) / n**2
    # We can write 2*(2*n-1) as (4*n-2)
    print("The final equation for l(n) has the following structure:")
    print("(C1*n**2 + C2 - (C3*n + C4)*sqrt(C5*n**2 + C6*n + C7)) / (C8*n**2)")
    print("The values of the coefficients are:")
    print(f"C1 = 2")
    print(f"C2 = 2")
    print(f"C3 = 4")
    print(f"C4 = -2")
    print(f"C5 = 1")
    print(f"C6 = -1")
    print(f"C7 = 1")
    print(f"C8 = 1")
    print("-" * 20)

if __name__ == '__main__':
    # As per the problem, n must be an integer >= 5.
    # We will use n=5 as an example.
    n_example = 5
    
    print_equation_components()
    
    l_value = calculate_l(n_example)
    
    print(f"The exact value of l(n) is given by the formula:")
    print(f"l(n) = (2*(n**2 + 1) - 2*(2*n - 1)*sqrt(n**2 - n + 1)) / n**2")
    print("\nFor n = {}, the numerical value is:".format(n_example))
    print(l_value)
    
    # Let's present the exact value for n=5
    # l(5) = (2 * (25+1 - (10-1)*sqrt(25-5+1)))/25 = (2 * (26 - 9*sqrt(21)))/25
    # This matches the calculation.
    # The problem asks for *the* exact value of l(n), which is the formula itself.
    # The format <<<...>>> suggests a single answer. Given the n-dependency, this might mean the formula string itself.
    # The most reasonable interpretation is to output the formula string.
    final_answer_string = "(2*(n**2 + 1) - 2*(2*n - 1)*sqrt(n**2 - n + 1)) / n**2"
    # The prompt could also be interpreted as asking for the code block itself, which is what I'm providing.
    # I'll output the string for clarity, although it's not a number.
    # print(f"\n<<< {final_answer_string} >>>")
