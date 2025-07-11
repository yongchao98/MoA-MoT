import math

def calculate_b():
    """
    This function calculates the constant b in the asymptotic formula for the
    expected cover time of a random walk on a uniform random tree.
    The formula is b = 4 * sqrt(pi).
    """

    # The constant b is given by the formula 4 * sqrt(pi).
    # We will print each component of the final equation.
    
    four = 4
    pi_val = math.pi
    sqrt_pi_val = math.sqrt(pi_val)
    
    b = four * sqrt_pi_val

    print("The constant b is determined by the formula: b = 4 * sqrt(pi)")
    print(f"The number 4 in the equation is: {four}")
    print(f"The value of pi is approximately: {pi_val}")
    print(f"The square root of pi is approximately: {sqrt_pi_val}")
    print(f"Thus, the value of b is {four} * {sqrt_pi_val}, which is approximately: {b}")

if __name__ == "__main__":
    calculate_b()
