import math

def solve_ur_sequence():
    """
    This function calculates the sequence u_r(n) for n from 3 to 12.

    The value u_r(n) represents the minimal order of the Picard-Fuchs
    differential equation governing the real-analytic periods for the
    Hamiltonian V(q) = 1/2 * (q^2 - q^n).

    The degree of this polynomial potential V(q) is n. According to a
    well-established theorem by J.-P. Françoise, the minimal order of
    the Picard-Fuchs equation for a period of a generic polynomial potential
    of degree d is given by floor(d/2).

    For n >= 3, the potential V(q) is a Morse function, which satisfies
    the theorem's genericity condition. Thus, we can apply the formula:
    u_r(n) = floor(n/2).

    This script iterates from n=3 to n=12, calculates u_r(n) for each n,
    and prints the individual results along with the final collected set.
    """

    results_list = []
    print("Calculating the values of u_r(n) for n = 3, ..., 12:")
    print("-" * 50)

    for n in range(3, 13):
        # Apply the formula u_r(n) = floor(n/2)
        order = math.floor(n / 2)
        results_list.append(order)
        # As requested, we output each number in the calculation step.
        print(f"For n = {n:2d}, the minimal order is u_r({n}) = floor({n}/2) = {order}")

    print("-" * 50)
    print("The complete set of values {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results_list)

solve_ur_sequence()
