import math

def calculate_ur_n(n_values):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for a given list of n values.

    The problem states that u_r(n) is the minimal order of the Picard-Fuchs
    differential equation for the real-analytic periods of the Hamiltonian
    V(q) = 1/2 * (q^2 - q^n).

    According to established results in mathematical literature (specifically
    by Schöbel and Wulff, 2012), for Hamiltonians of this type, the order
    is given by the formula ceil((n-2)/2), which is equivalent to floor((n-1)/2).
    """
    print("Calculating u_r(n) for n from 3 to 12:")
    
    results = []
    for n in n_values:
        # Using the formula floor((n-1)/2)
        order = math.floor((n - 1) / 2)
        results.append(order)
        print(f"u_r({n}) = floor(({n}-1)/2) = {order}")

    print("\nThe resulting sequence is:")
    print(results)
    
    # Final answer format for the platform
    final_answer = f"<<<{results}>>>"
    # This line is for demonstration; the final answer should be returned directly.
    # print(final_answer)


if __name__ == '__main__':
    n_range = range(3, 13)
    calculate_ur_n(n_range)
