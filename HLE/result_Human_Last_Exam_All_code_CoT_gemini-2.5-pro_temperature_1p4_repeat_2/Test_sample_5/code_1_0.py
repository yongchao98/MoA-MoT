import sys

def solve_proportionality_factor():
    """
    Calculates the proportionality factor in the gamma matrix identity.

    The identity is:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}

    The user is prompted to enter the dimension 'd' and the rank 'k'.
    """
    try:
        d_str = input("Enter the dimension of spacetime, d: ")
        d = int(d_str)
        if d < 1:
            print("Dimension d must be a positive integer.")
            return

        k_str = input("Enter the rank of the gamma matrix, k: ")
        k = int(k_str)
        if not (0 <= k <= d):
            print(f"Rank k must be an integer between 0 and d={d}.")
            return

    except ValueError:
        print("Invalid input. Please enter integers for d and k.")
        return

    # The derived proportionality factor C(d, k) = 4k^2 - 4dk + d - d^2
    factor = 4 * k**2 - 4 * d * k + d - d**2

    print("\nThe relation is: gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}")
    print(f"\nFor d = {d} and k = {k}, the proportionality factor C(d, k) is derived from the formula:")
    print("C(d, k) = 4*k^2 - 4*d*k + d - d^2\n")
    print(f"Plugging in the values:\nC({d}, {k}) = 4*({k})^2 - 4*({d})*({k}) + ({d}) - ({d})^2")
    
    term1 = 4 * k**2
    term2 = -4 * d * k
    term3 = d
    term4 = -d**2
    
    # To avoid printing -- for negative numbers
    if term2 >= 0:
      s_term2 = f"+ {term2}"
    else:
      s_term2 = f"- {-term2}"
    if term3 >= 0:
      s_term3 = f"+ {term3}"
    else:
      s_term3 = f"- {-term3}"
    if term4 >= 0:
      s_term4 = f"+ {term4}"
    else:
      s_term4 = f"- {-term4}"

    print(f"         = {term1} {s_term2} {s_term3} {s_term4}")
    print(f"         = {factor}")

    # The final answer in the required format
    print(f"\n<<<The proportionality factor is {factor}>>>")

if __name__ == '__main__':
    solve_proportionality_factor()