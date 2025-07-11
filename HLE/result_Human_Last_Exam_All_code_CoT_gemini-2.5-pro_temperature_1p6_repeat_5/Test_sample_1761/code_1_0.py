import sys
import math

def main():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensor O_{P^n}(2) on complex projective space P^n.
    """
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <n>")
        print("where <n> is the dimension of the projective space.")
        return

    try:
        n = int(sys.argv[1])
        if n < 1:
            raise ValueError
    except ValueError:
        print("Error: Please provide a positive integer for n.")
        return

    print(f"Calculating the dimension for n = {n}:")
    
    # h^0(P^n, O(k)) = C(n+k, k)
    # math.comb(n, k) computes the binomial coefficient "n choose k"

    # h^0(P^n, O(1)^{n+1}) = (n+1) * h^0(P^n, O(1))
    h0_O1 = math.comb(n + 1, 1)
    h0_O1_oplus = (n + 1) * h0_O1
    
    # h^0(P^n, O(2))
    h0_O2 = math.comb(n + 2, 2)
    
    # The dimension is the difference: h^0(Omega^1(2)) = h^0(O(1)^{n+1}) - h^0(O(2))
    result = h0_O1_oplus - h0_O2
    
    print("\nBased on the Euler sequence, the dimension is calculated as:")
    print(f"h^0(Omega^1(2)) = h^0(O(1)^{n+1}) - h^0(O(2))")
    
    print("\nFirst, we compute the dimensions of the known spaces:")
    print(f"h^0(O(1)^{n+1}) = {h0_O1_oplus}")
    print(f"h^0(O(2)) = {h0_O2}")
    
    print("\nSubstituting these values, we get the final equation:")
    print(f"{h0_O1_oplus} - {h0_O2} = {result}")

    print(f"\nThe complex dimension is {result}.")
    print(f"Note that this result is equal to C(n+1, 2) = C({n+1}, 2).")


if __name__ == "__main__":
    main()
