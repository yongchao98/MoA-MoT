import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_{P^n} tensor O_{P^n}(2).
    """
    # Please set the value of n for the projective space P^n.
    # n must be a non-negative integer. Example is set to n=3.
    n = 3

    print(f"Calculating the dimension for the complex projective space P^n where n = {n}.")

    # The dimension is given by the formula:
    # h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^{n+1}) - h^0(P^n, O(2))

    # Step 1: Calculate h^0(P^n, O(1)^{n+1}) = (n+1) * C(n+1, 1) = (n+1)^2
    h0_O1_sum = (n + 1)**2

    # Step 2: Calculate h^0(P^n, O(2)) = C(n+2, 2)
    h0_O2 = math.comb(n + 2, 2)

    # Step 3: Calculate the final dimension
    dimension = h0_O1_sum - h0_O2

    print("\nCalculation steps:")
    print(f"h^0(P^{n}, O(1)^{{n+1}}) = ({n}+1)^2 = {h0_O1_sum}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")
    print(f"The dimension of H^0(P^{n}, Omega^1(2)) is the difference:")
    print(f"Dimension = {h0_O1_sum} - {h0_O2} = {dimension}")

    # The result can also be calculated directly with the formula C(n+1, 2)
    direct_formula_result = math.comb(n + 1, 2)
    print(f"\nThis confirms the simplified formula: C(n+1, 2) = {direct_formula_result}")

if __name__ == '__main__':
    calculate_dimension()