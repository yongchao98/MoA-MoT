import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n(2) on complex projective space P^n.
    """
    try:
        n_str = input("Please enter the dimension n of the complex projective space P^n_C (n >= 1): ")
        n = int(n_str)
        if n < 1:
            print("The dimension n must be an integer greater than or equal to 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Calculate the dimension of H^0(P^n, O(1)^{n+1})
    # h^0(P^n, O(k)) = binomial(n+k, k)
    # h^0(P^n, O(1)) = binomial(n+1, 1) = n+1
    h0_O1_bundle_val = (n + 1) * (n + 1)

    # Calculate the dimension of H^0(P^n, O(2))
    # h^0(P^n, O(2)) = binomial(n+2, 2)
    h0_O2_val = math.comb(n + 2, 2)

    # The result is the difference, based on the short exact sequence of global sections
    # derived from the twisted Euler sequence.
    result = h0_O1_bundle_val - h0_O2_val

    print("\nTo find the dimension of the space of global sections, we use the formula:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^{n+1}) - h^0(P^n, O(2))")
    print("\nFor n = {}, we calculate the terms:".format(n))
    print("h^0(P^{}, O(1)^{{{}+1}}) = ({}+1)^2 = {}".format(n, n, n, h0_O1_bundle_val))
    print("h^0(P^{}, O(2)) = C({}+2, 2) = {}".format(n, n, h0_O2_val))
    print("\nThe dimension is the difference:")
    print("dim = {} - {} = {}".format(h0_O1_bundle_val, h0_O2_val, result))
    print("\nThis simplifies to the formula n(n+1)/2, which for n={} is {}.".format(n, math.comb(n+1,2)))


if __name__ == '__main__':
    calculate_dimension()
