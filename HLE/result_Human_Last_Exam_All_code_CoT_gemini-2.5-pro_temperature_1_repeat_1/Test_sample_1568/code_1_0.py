import numpy as np
from scipy.special import gamma

def calculate_infinite_product(z):
    """
    Calculates the value of the infinite product Prod_{n=3 to inf} (1 - z^3/n^3).

    The function uses a closed-form expression involving the Gamma function:
    1 / ( (1 - z^3/1) * (1 - z^3/8) * Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2) )
    where w is the complex cube root of unity.

    Args:
        z (complex or float): The value of z in the expression.

    Returns:
        complex: The result of the infinite product.
    """
    # Define the principal cube root of unity, omega = e^(i*2*pi/3)
    omega = np.exp(2j * np.pi / 3)
    # The other non-real root is omega^2
    omega2 = omega**2

    # The product starts from n=3. We need to divide by the terms for n=1 and n=2.
    n_start = 3
    
    # Calculate the term for n=1: (1 - z^3/1^3)
    term_n1_val = 1**3
    term1 = 1 - z**3 / term_n1_val
    
    # Calculate the term for n=2: (1 - z^3/2^3)
    term_n2_val = 2**3
    term2 = 1 - z**3 / term_n2_val

    # Calculate the product of the Gamma functions based on the known identity.
    gamma_prod = gamma(1 - z) * gamma(1 - z * omega) * gamma(1 - z * omega2)

    # The full infinite product from n=1 is 1 / gamma_prod.
    # We divide by the first two terms to get the product from n=3.
    result = 1 / (term1 * term2 * gamma_prod)

    # Output the components of the final equation
    print(f"Calculating the product for z = {z}")
    print(f"The infinite product starts from n = {n_start}")
    print("The formula is: 1 / [ (1 - z^3/{}) * (1 - z^3/{}) * G(1-z) * G(1-z*w) * G(1-z*w^2) ]".format(term_n1_val, term_n2_val))
    print("-" * 30)
    print(f"Term for n=1: (1 - {z}^3 / {term_n1_val}) = {term1}")
    print(f"Term for n=2: (1 - {z}^3 / {term_n2_val}) = {term2}")
    print(f"Product of Gamma functions part: G(1-{z})*G(1-{z}*w)*G(1-{z}*w^2) = {gamma_prod}")
    print("-" * 30)
    print(f"Final Result: {result}")
    
    return result

if __name__ == '__main__':
    # Example: Calculate the product for z = 0.5
    # The value of z must not be an integer, or an integer multiple of omega or omega^2.
    z_value = 0.5
    calculate_infinite_product(z_value)