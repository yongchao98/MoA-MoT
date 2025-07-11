import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_P^n(2) on complex projective n-space.
    """
    try:
        n_str = input("Enter the dimension n of the projective space P^n: ")
        n = int(n_str)
        if n < 1:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The dimension h^0(O(k)) is given by the binomial coefficient C(n+k, k).
    # We use the relationship: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2)).
    
    # Calculate h^0(O(1))
    h0_O1 = math.comb(n + 1, 1)
    
    # Calculate h^0(O(1)^(n+1)) = (n+1) * h^0(O(1))
    h0_O1_bundle = (n + 1) * h0_O1

    # Calculate h^0(O(2))
    h0_O2 = math.comb(n + 2, 2)

    # The final dimension is the difference
    result = h0_O1_bundle - h0_O2
    
    # Print the step-by-step calculation
    print("\n--- Calculation Breakdown ---")
    print(f"For n = {n}:")
    print(f"The dimension h^0(P^{n}, O(1)) is calculated as C({n}+{1}, 1) = {h0_O1}.")
    print(f"The dimension h^0(P^{n}, O(1)^({n}+1)) is therefore ({n}+1) * {h0_O1} = {h0_O1_bundle}.")
    print(f"The dimension h^0(P^{n}, O(2)) is calculated as C({n}+{2}, 2) = {h0_O2}.")
    print("\nThe final equation for the dimension is: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
    print(f"Result: {h0_O1_bundle} - {h0_O2} = {result}")
    
    # Verify with the simplified formula n(n+1)/2
    verification = n * (n + 1) // 2
    print(f"\nThis matches the simplified general formula n*(n+1)/2 = {verification}.")


if __name__ == "__main__":
    calculate_dimension()
