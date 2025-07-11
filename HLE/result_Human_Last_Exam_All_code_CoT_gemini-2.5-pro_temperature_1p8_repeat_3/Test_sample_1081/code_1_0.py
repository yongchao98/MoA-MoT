import math

def calculate_e8_tori_number():
    """
    Calculates the number of Fq-rational maximal tori of a reductive group of type E_8.
    This number is equal to the order of the Weyl group W(E_8).
    """

    # The order of the Weyl group W(E_8) is given by the product:
    # 2^14 * 3^5 * 5^2 * 7
    p1 = 2
    e1 = 14
    n1 = p1**e1

    p2 = 3
    e2 = 5
    n2 = p2**e2

    p3 = 5
    e3 = 2
    n3 = p3**e3

    p4 = 7
    e4 = 1
    n4 = p4**e4

    # Calculate the total number
    total_number = n1 * n2 * n3 * n4

    # Print the explanation and the calculation
    print("The number of F_q-rational maximal tori in a group of type E_8 is given by the order of its Weyl group, |W(E_8)|.")
    print("The formula for the order is: 2^14 * 3^5 * 5^2 * 7.")
    print("\nLet's calculate the value of each term in the final equation:")
    print(f"{p1}^{e1} = {n1}")
    print(f"{p2}^{e2} = {n2}")
    print(f"{p3}^{e3} = {n3}")
    print(f"{p4}^{e4} = {n4}")

    print("\nNow, we multiply these numbers together to get the final result.")
    print(f"The final equation is: {n1} * {n2} * {n3} * {n4} = {total_number}")

    print(f"\nThus, the exact number of F_q-rational maximal tori of G is {total_number}.")

# Execute the function
calculate_e8_tori_number()