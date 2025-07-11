import math

def calculate_grassmann_integral(a, b):
    """
    This function demonstrates the calculation of a simple Grassmann integral:
    Z = ∫ dη exp(a + bη)
    based on the fundamental rules of Berezin integration.

    The Pauli exclusion principle is encoded in the property η^2 = 0.
    This means the Taylor series for exp(bη) is simply 1 + bη.

    The Berezin integral is defined by two rules:
    1. ∫ dη = 0
    2. ∫ η dη = 1
    """

    print("Task: Calculate the Grassmann integral Z = ∫ dη exp(a + bη)")
    print(f"Given a = {a} and b = {b}\n")

    print("Step 1: Expand the exponential using η^2 = 0.")
    print("The Pauli exclusion principle implies that for a Grassmann variable η, η^2 = 0.")
    print("Therefore, the Taylor expansion of the exponential term exp(bη) truncates after the linear term:")
    print("exp(bη) = 1 + bη + (bη)^2/2! + ... = 1 + bη")
    print(f"So, the full expression becomes:")
    print(f"exp({a} + {b}η) = exp({a}) * exp({b}η) = exp({a}) * (1 + {b}η)")
    exp_a = math.exp(a)
    print(f"where exp({a}) ≈ {exp_a:.4f}\n")


    print("Step 2: Apply the rules of Berezin integration to the expanded form.")
    print("Z = ∫ dη [exp(a) * (1 + bη)]")
    print("By linearity, we can write this as:")
    print("Z = exp(a) * [∫ dη * 1 + b * ∫ dη * η]")
    print("Now, we apply the definitional rules of the integral:")
    # Rule 1: ∫ dη * 1 = 0
    integral_of_1 = 0
    # Rule 2: ∫ dη * η = 1
    integral_of_eta = 1
    print(f"Rule 1: ∫ dη = {integral_of_1}")
    print(f"Rule 2: ∫ η dη = {integral_of_eta}")
    print(f"Substituting these values into the equation for Z:")
    print(f"Z = exp({a}) * [{integral_of_1} + {b} * {integral_of_eta}]")

    print("\nStep 3: Calculate the final numerical result.")
    result = exp_a * (integral_of_1 + b * integral_of_eta)
    print(f"Z = {exp_a:.4f} * ({integral_of_1} + {b})")
    print(f"Z = {exp_a:.4f} * {b}")
    print(f"Z = {result:.4f}")
    
    print("\n--- Final Equation ---")
    print("The symbolic result is: ∫ dη exp(a + bη) = b * exp(a)")
    print(f"Substituting the given numbers: ∫ dη exp({a} + {b}η) = {b} * {math.exp(a):.4f} = {result:.4f}")


# Example usage with arbitrary numbers for the coefficients a and b
calculate_grassmann_integral(a=2.0, b=5.0)