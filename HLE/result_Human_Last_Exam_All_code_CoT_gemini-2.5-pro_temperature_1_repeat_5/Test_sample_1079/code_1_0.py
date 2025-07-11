import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group that have a regular
    eigenvector with an eigenvalue of order 10.
    """
    
    # Parameters for the reflection group of type H3
    group_order = 120  # |W|
    coxeter_number = 10 # h

    # An element has a regular eigenvector iff it's a regular element.
    # A regular element in H3 has an eigenvalue of order h=10.
    # So, we need to find the number of regular elements in H3.
    # The formula is N = |W| * phi(h) / h, where phi is Euler's totient function.

    # Calculate phi(10)
    # Integers coprime to 10 are 1, 3, 7, 9. There are 4 such numbers.
    phi_h = 0
    for i in range(1, coxeter_number + 1):
        if math.gcd(i, coxeter_number) == 1:
            phi_h += 1

    # Calculate the number of regular elements
    num_regular_elements = group_order * phi_h // coxeter_number

    # Print the final equation as requested
    print("The number of such elements is calculated by the formula for regular elements in a Coxeter group:")
    print("Number = |W| * phi(h) / h")
    print(f"For H3, the group order |W| is {group_order} and the Coxeter number h is {coxeter_number}.")
    print(f"Euler's totient function phi({coxeter_number}) is {phi_h}.")
    print("The final calculation is:")
    print(f"{group_order} * {phi_h} / {coxeter_number} = {num_regular_elements}")

solve_h3_regular_elements()