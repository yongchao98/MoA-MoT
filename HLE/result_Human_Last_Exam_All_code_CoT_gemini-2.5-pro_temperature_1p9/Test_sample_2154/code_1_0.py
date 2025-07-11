import math

def calculate_ur(n_list):
    """
    Calculates the minimal order u_r(n) for a list of n values.
    The formula is u_r(n) = floor((n-1)/2), based on the theory of
    Picard-Fuchs equations for periods of Hamiltonian systems with
    a polynomial potential of degree n.
    """
    results = {}
    for n in n_list:
        # The potential V(q) = 1/2 * (q^2 - q^n) has degree n.
        # The minimal order of the Picard-Fuchs equation is floor((n-1)/2).
        order = math.floor((n - 1) / 2)
        results[n] = order
    return results

def main():
    """
    Main function to calculate and print the sequence {u_r(3), u_r(4), ..., u_r(12)}.
    """
    n_values = range(3, 13)
    ur_values = calculate_ur(n_values)

    print("The required sequence is {u_r(3), u_r(4), ..., u_r(12)}.")
    print("The minimal order u_r(n) is calculated using the formula: floor((n-1)/2).")
    print("-" * 30)

    # Print each calculation as an "equation"
    output_list = []
    for n, val in ur_values.items():
        print(f"For n = {n}, u_r({n}) = floor(({n}-1)/2) = {val}")
        output_list.append(str(val))

    print("-" * 30)
    print("The complete set is:")
    print("{" + ", ".join(output_list) + "}")


if __name__ == "__main__":
    main()