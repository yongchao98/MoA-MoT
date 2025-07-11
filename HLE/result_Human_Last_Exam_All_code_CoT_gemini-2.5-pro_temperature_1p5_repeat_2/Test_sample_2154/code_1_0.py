def calculate_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n).

    Args:
      n: An integer representing the power in the potential V(q).

    Returns:
      The minimal order u_r(n).
    """
    if n % 2 == 1:
        # For odd n, the order is n - 1
        return n - 1
    else:
        # For even n, the order is n / 2
        return n // 2

def main():
    """
    Calculates and prints u_r(n) for n from 3 to 12.
    """
    print("The minimal orders u_r(n) for n = 3 to 12 are:")
    results = []
    for n in range(3, 13):
        order = calculate_ur(n)
        results.append(order)
        print(f"u_r({n}) = {order}")
    
    # The problem asks for the set of values. Let's print it in a set-like format.
    # Note that the order corresponds to n=3, 4, ..., 12.
    result_str = ", ".join(map(str, results))
    print(f"\nThe sequence {{u_r(3), ..., u_r(12)}} is: {{{result_str}}}")

if __name__ == "__main__":
    main()