def print_formula():
    """
    Prints the derived exact formula for l_k(n).
    """
    # The derived formula for l_k(n) is:
    # l_k(n) = (1/2) * ln(n + 1) - (2 - 1/n) * k^2 + (n - 1) * ln(k)
    # The code below prints this formula as a string.
    
    formula_string = "l_k(n) = (1/2) * ln(n + 1) - (2 - 1/n) * k^2 + (n - 1) * ln(k)"
    
    print("The exact value of l_k(n) in terms of n and k is:")
    print(formula_string)

if __name__ == '__main__':
    print_formula()