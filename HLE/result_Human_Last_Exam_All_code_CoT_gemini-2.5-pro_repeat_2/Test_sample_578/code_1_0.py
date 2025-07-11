import math

def get_d(n):
    """
    Calculates the dimension d_n based on the provided formula.
    d_1 = 2
    d_n = C_{n-1} for n > 1, where C_k is the k-th Catalan number.
    """
    if n == 1:
        return 2
    elif n > 1:
        # Calculate the (n-1)-th Catalan number C_{n-1}
        k = n - 1
        # C_k = (1/(k+1)) * comb(2k, k)
        return math.comb(2 * k, k) // (k + 1)
    else:
        raise ValueError("n must be a positive integer")

def solve():
    """
    Calculates and prints the product of d_n for n from 1 to 8.
    """
    d_values = []
    print("Calculating the values of d_n for n = 1 to 8:")
    for n in range(1, 9):
        d_val = get_d(n)
        d_values.append(d_val)
        print(f"d_{n} = {d_val}")

    product = 1
    for val in d_values:
        product *= val
    
    equation_str = " * ".join(map(str, d_values))
    print("\nThe final product is calculated as:")
    print(f"d_1 * d_2 * d_3 * d_4 * d_5 * d_6 * d_7 * d_8 = {equation_str} = {product}")

if __name__ == "__main__":
    solve()