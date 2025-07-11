import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_{P^n} tensored by O_{P^n}(2).

    The dimension is given by the formula:
    h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))
    """
    try:
        n_str = input("Enter the dimension n of the complex projective space P^n (n>=1): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Calculate h^0(P^n, O(1)^(n+1))
    # h^0(P^n, O(1)) = C(n+1, 1) = n+1
    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
    h0_O1_oplus = (n + 1) * (n + 1)
    
    # Calculate h^0(P^n, O(2))
    # h^0(P^n, O(2)) = C(n+2, 2)
    h0_O2 = math.comb(n + 2, 2)

    # The result is the difference
    result = h0_O1_oplus - h0_O2

    print("\nStarting the calculation for n = {}...".format(n))
    print("The dimension is the kernel of a map in the long exact sequence of cohomology.")
    print("dim = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print("-" * 30)
    print("First term: h^0(P^n, O(1)^(n+1)) = (n+1) * C(n+1, 1)")
    print("h^0(P^{}, O(1)^(n+1)) = ({} + 1) * C({} + 1, 1) = {} * {} = {}".format(n, n, n, n + 1, n + 1, h0_O1_oplus))
    print("-" * 30)
    print("Second term: h^0(P^n, O(2)) = C(n+2, 2)")
    print("h^0(P^{}, O(2)) = C({} + 2, 2) = C({}, 2) = {}".format(n, n, n + 2, h0_O2))
    print("-" * 30)
    print("Final equation:")
    print("{} - {} = {}".format(h0_O1_oplus, h0_O2, result))
    print("-" * 30)
    print("The complex dimension is {}.".format(result))

calculate_dimension()