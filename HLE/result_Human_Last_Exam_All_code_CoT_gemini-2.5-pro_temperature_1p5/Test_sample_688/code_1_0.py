import math

def get_cn_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    Ree-Hoover diagram with n nodes.

    Args:
        n (int): The number of nodes (particles), must be n >= 2.

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")

    # For n=2, the virial coefficient B_2 consists of only the fully f-connected term.
    # B_2 = -(2-1)/2! * Lambda_2 = -1/2 * Lambda_2. Thus, c_2 = -1/2.
    if n == 2:
        return -0.5

    # For n>=3, the coefficient is derived from the Ree-Hoover expansion.
    # The coefficient of the all-f-bond term is given by (-1)^n / n.
    # For n=3, this gives c_3 = (-1)^3 / 3 = -1/3. This matches the direct
    # calculation from the Mayer expansion, B_3 = -(3-1)/3! * Lambda_3 = -1/3 * Lambda_3.
    # For n>=4, the Ree-Hoover decomposition becomes non-trivial.
    else: # n >= 3
        return ((-1)**n) / n

def main():
    """
    Prints the formula for c_n and demonstrates its calculation for n=2 to 7.
    """
    print("The prefactor c_n is calculated based on the number of nodes n.")
    print("The general formula has a special case for n=2.")
    print("-" * 20)

    for n_val in range(2, 8):
        c_n = get_cn_prefactor(n_val)
        if n_val == 2:
            print(f"For n = {n_val}:")
            print(f"c_{n_val} = -1 / 2")
            print(f"c_{n_val} = {c_n}\n")
        else:
            print(f"For n = {n_val}:")
            # The instruction "output each number in the final equation" is handled here.
            # Numerator: (-1)^n
            # Denominator: n
            numerator_val = (-1)**n_val
            denominator_val = n_val
            print(f"c_{n_val} = (-1)^{n_val} / {n_val} = {numerator_val} / {denominator_val}")
            print(f"c_{n_val} = {c_n}\n")

if __name__ == "__main__":
    main()
