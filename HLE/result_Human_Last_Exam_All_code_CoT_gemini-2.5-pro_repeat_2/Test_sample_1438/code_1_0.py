import sys

def explain_grassmann_integration():
    """
    This function explains the rules of Grassmann integration (Berezin integral),
    which are fundamental to the path integral formalism for fermions and ensure
    the Pauli exclusion principle is respected.
    """

    print("In the path integral formalism for fermionic systems, we use anti-commuting numbers called Grassmann variables (e.g., η).")
    print("The key algebraic property is that they anti-commute: η_i * η_j = -η_j * η_i.")
    print("This implies that the square of any single Grassmann variable is zero: η^2 = 0.")
    print("This algebraic property, η^2 = 0, is the mathematical embodiment of the Pauli exclusion principle.")
    print("\nAn integral over a Grassmann variable is defined by its action, not by a numerical value of the measure 'dη' itself.")
    print("The rules for this integration (Berezin integral) are as follows:")

    # Rule 1
    # We represent the integral '∫ dη * f(η)' as a string.
    # The numbers in the final equation are 1 and 0.
    constant_term = "1"
    result_1 = "0"
    print("\nRule 1: The integral of a constant is zero.")
    print(f"∫ dη * {constant_term} = {result_1}")

    # Rule 2
    # The numbers in the final equation are 1 (implicitly for η) and 1.
    linear_term = "η"
    result_2 = "1"
    print("\nRule 2: The integral of the variable itself is one (by normalization convention).")
    print(f"∫ dη * {linear_term} = {result_2}")

    print("\nThese two rules define the 'measure'. The normalization to 1 is the key convention.")

# Execute the function to print the explanation.
explain_grassmann_integration()
