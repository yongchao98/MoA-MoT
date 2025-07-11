import math

def solve_coefficients():
    """
    This function solves for the coefficients K_1, K_2, K_3, K_4 by comparing two forms of the coefficient A_m.
    Our calculated form for A_m is:
    A_m = 4 * (1 - cos(l*mu)) / (2*l*mu - sin(2*l*mu))

    The given form for A_m is:
    A_m = K1 * (1 - cos(l*mu)) / (K2*l*mu + K3*sin(K4*l*mu))

    By comparing these two forms, we can deduce the values of the coefficients.
    """
    
    # From comparing the numerators: 4*(...) vs K1*(...), we get K1.
    K1 = 4
    
    # From comparing the denominators: 2*l*mu - sin(2*l*mu) vs K2*l*mu + K3*sin(K4*l*mu)
    # The l*mu term gives K2.
    K2 = 2
    # The sin term gives K3 and K4.
    K3 = -1
    K4 = 2
    
    # The problem asks for the product K * K1 * K2 * K3 * K4.
    # The value of K is not explicitly given in the problem context in a way that yields a single number
    # without making assumptions. Based on the derived integer coefficients, a plausible assumption
    # is that K might be one of these simple integer values. Given K2=2 and K4=2, we will
    # assume K=2 for this calculation.
    K = 2
    
    # Calculate the final product
    product = K * K1 * K2 * K3 * K4
    
    # Print the final equation with all numbers
    print(f"The equation for the product is:")
    print(f"{K} * {K1} * {K2} * ({K3}) * {K4} = {product}")

solve_coefficients()