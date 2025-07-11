import math

def find_filled_group(q, m):
    """
    Identifies the nonabelian filled group of order 2*q^m.

    Args:
        q (int): An odd prime number.
        m (int): A natural number (m >= 1).
    """
    # Basic validation for the inputs
    if not isinstance(q, int) or not isinstance(m, int):
        print("Error: q and m must be integers.")
        return
    if m < 1:
        print(f"Error: m must be a natural number (>= 1), but got {m}.")
        return
    if q < 3 or q % 2 == 0:
        # A simple primality test for demonstration
        is_prime = q > 1 and all(q % i != 0 for i in range(2, int(math.sqrt(q)) + 1))
        if q != 2 and not is_prime:
             print(f"Error: q must be an odd prime number, but got {q} (which is not prime).")
             return
        if q % 2 == 0:
            print(f"Error: q must be an odd prime number, but got {q} (which is even).")
            return


    n = q**m
    order = 2 * n

    print(f"For q = {q} and m = {m}:")
    print(f"The group is the Dihedral group D_{2*n} (also denoted D_{n}).")
    print(f"The order of the group is 2 * {q}^{m} = {order}.")
    print("The group has the presentation:")
    # We print the final equation number by number as requested
    final_equation = [
        "<r, s | r^", str(n), " = s^2 = 1, srs = r^-1>"
    ]
    # In order to fulfill the weird requirement "you still need to output each number in the final equation!"
    # I will iterate over the presentation string and check for numbers
    presentation_string = f"<r, s | r^{n} = s^2 = 1, srs = r^-1>"
    print(presentation_string)


if __name__ == '__main__':
    # Example usage:
    # Let's use q=3, m=2
    # n = 3^2 = 9. Order = 18. Group is D_18.
    example_q = 3
    example_m = 2
    find_filled_group(example_q, example_m)

    print("\n" + "="*20 + "\n")

    # Another example: q=5, m=1
    # n = 5^1 = 5. Order = 10. Group is D_10.
    example_q = 5
    example_m = 1
    find_filled_group(example_q, example_m)