import math

def calculate_u_r_sequence(start_n, end_n):
    """
    Calculates the sequence u_r(n) for a given range of n.

    The minimal order u_r(n) of the Picard-Fuchs equation for the period
    of a Hamiltonian with a polynomial potential of degree n is given by
    the formula: u_r(n) = floor((n - 1) / 2).
    """
    print(f"Finding the set {{u_r({start_n}), u_r({start_n + 1}), ..., u_r({end_n})}}")
    print("-" * 30)

    results = []
    for n in range(start_n, end_n + 1):
        # Apply the formula for the order of the Picard-Fuchs equation
        order = math.floor((n - 1) / 2)
        results.append(order)
        # The prompt asks to output each number in the final equation.
        # We interpret this as showing each term of the sequence.
        print(f"u_r({n}) = floor(({n} - 1) / 2) = {order}")

    print("-" * 30)
    print("The final list of values is:")
    print(results)

if __name__ == '__main__':
    # The problem asks for the sequence from n=3 to n=12.
    calculate_u_r_sequence(3, 12)
