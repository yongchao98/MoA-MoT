def print_asymptotic_formula():
    """
    Prints the asymptotic formula for d_{B,delta}.
    The formula uses B and L, where L = log(delta^{-1}).
    The formula is a combination of the two main asymptotic regimes.
    """
    formula = "sqrt(B*L) + L/log(L)"
    print(formula)

if __name__ == '__main__':
    print_asymptotic_formula()