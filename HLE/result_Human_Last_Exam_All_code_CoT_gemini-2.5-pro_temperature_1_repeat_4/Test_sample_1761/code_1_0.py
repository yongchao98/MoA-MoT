import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on n-dimensional complex projective space P^n.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n: ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension 'n' must be a positive integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for 'n'.")
        return

    print(f"\nCalculating the dimension for n = {n}:")
    print("The dimension is given by the formula:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print("h^0(P^n, Omega^1(2)) = (n+1)^2 - C(n+2, 2)\n")

    # Step 1: Calculate (n+1)^2
    term1_base = n + 1
    term1 = term1_base**2
    
    # Step 2: Calculate C(n+2, 2)
    comb_n = n + 2
    comb_k = 2
    term2 = math.comb(comb_n, comb_k)

    # Step 3: Calculate the final result
    result = term1 - term2

    # Print the detailed calculation
    print("Plugging in n = {}:".format(n))
    print(f"dim = ({n} + 1)^2 - C({n} + 2, 2)")
    print(f"dim = {term1_base}^2 - C({comb_n}, {comb_k})")
    print(f"dim = {term1} - {term2}")
    print(f"dim = {result}")

    # Also, we can show the simplified formula result
    simplified_result = n * (n + 1) // 2
    print(f"\nUsing the simplified formula C(n+1, 2) = C({n+1}, 2) = {simplified_result}")


if __name__ == "__main__":
    calculate_dimension()