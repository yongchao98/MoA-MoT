import math

def solve_expected_T():
    """
    This function calculates and prints the expected value of T.
    The derivation shows that E[T] = 3 - 2 * e^(1/3).
    """

    # The final equation is E[T] = A - B * e^(C)
    A = 3
    B = 2
    C = 1/3

    # Calculate the numerical value
    e_to_C = math.exp(C)
    expected_T = A - B * e_to_C

    # Print the explanation and results
    print("The expected value E[T] is derived to be of the form: A - B * e^(C)")
    print("The numbers in this final equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = 1/{A}") # A is 3, so this is 1/3

    print(f"\nThe symbolic equation is: E[T] = {A} - {B} * e^(1/3)")
    print(f"The numerical value of e^(1/3) is: {e_to_C}")
    print(f"The final numerical value for E[T] is: {expected_T}")

# Execute the function to find the answer
solve_expected_T()