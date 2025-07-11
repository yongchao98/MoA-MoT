import sys

def solve_toric_code_degeneracy():
    """
    Calculates the ground space degeneracy of the toric code on a plane
    with a given number of smooth and rough holes.
    """
    # Example values for n (smooth holes) and m (rough holes)
    # You can change these values to test other scenarios.
    n = 3
    m = 2

    # The problem asks for the general formula. Based on the analysis, the
    # number of logical qubits k on an infinite plane is n + m.
    # The ground space degeneracy is 2^k.

    print(f"For a toric code on a plane with n={n} smooth holes and m={m} rough holes:")
    print(f"The number of logical qubits is k = n + m = {n} + {m} = {n+m}.")
    print("The ground space degeneracy is 2^k.")
    
    # The final answer is the formula itself.
    # We output the components of the formula as requested.
    print(f"Final Equation: 2^({n} + {m})")

solve_toric_code_degeneracy()