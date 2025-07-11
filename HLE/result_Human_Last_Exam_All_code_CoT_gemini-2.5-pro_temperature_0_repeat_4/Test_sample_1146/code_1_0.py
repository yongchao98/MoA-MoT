import math

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def find_n_for_ring_graph():
    """
    Finds the values of n for which AG(Z_n) is a "ring graph".
    Based on the analysis, this corresponds to finding n such that phi(n) = 2.
    """
    result_n = []
    # We only need to check a small range as phi(n) grows reasonably fast.
    for n in range(2, 50):
        if phi(n) == 2:
            result_n.append(n)
    return result_n

if __name__ == "__main__":
    solutions = find_n_for_ring_graph()
    # The problem asks to format the output as a sequence in a set.
    print(f"n in {{ {', '.join(map(str, solutions))} }}")
    # The final answer format for the platform
    print(f"<<<n in {{ {', '.join(map(str, solutions))} }}>>>")
