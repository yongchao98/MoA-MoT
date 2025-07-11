def berezin_integral(A, B):
    """
    Performs the Berezin integral for a function of a single Grassmann variable.
    A function of a single Grassmann variable eta can be written as f(eta) = A + B*eta.
    The integral is defined as: integral d(eta) f(eta) = B.

    Args:
        A (float or int): The constant part of the function.
        B (float or int): The coefficient of the eta term.

    Returns:
        float or int: The result of the integration, which is B.
    """
    return B

def main():
    """
    Demonstrates the rules of Berezin integration for two primary cases.
    """
    print("In Grassmann algebra, any analytic function of a variable 'eta' is f(eta) = A + B*eta, since eta*eta = 0.")
    print("The Berezin integral is defined by the rule: integral d(eta) * f(eta) = B.\n")

    # Case 1: Integrate a constant function, f(eta) = 1.
    # Here, A=1, B=0.
    A1, B1 = 1, 0
    result1 = berezin_integral(A1, B1)
    print(f"Case 1: Integrating f(eta) = {A1} + {B1} * eta")
    # Using 'η' for a more standard look in the output string.
    print(f"∫ d(η) ({A1} + {B1}η) = {result1}")
    print("-" * 20)

    # Case 2: Integrate the Grassmann variable itself, f(eta) = eta.
    # Here, A=0, B=1.
    A2, B2 = 0, 1
    result2 = berezin_integral(A2, B2)
    print(f"Case 2: Integrating f(eta) = {A2} + {B2} * eta")
    print(f"∫ d(η) ({A2} + {B2}η) = {result2}")
    print("-" * 20)
    
    print("\nThe second case, ∫d(η)η = 1, defines the normalization of the integration measure.")


if __name__ == "__main__":
    main()
