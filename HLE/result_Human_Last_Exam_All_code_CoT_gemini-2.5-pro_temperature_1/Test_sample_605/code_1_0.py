import sys

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau hypersurface.
    """
    # The weights of the ambient space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]

    # For a Calabi-Yau hypersurface, the degree 'd' of the defining polynomial
    # is equal to the sum of the weights.
    # We can calculate this sum.
    sum_of_weights = sum(weights)
    d = sum_of_weights

    # The Crawley-Nordström invariant is defined as c_N = d - sum(w_i).
    c_N = d - sum_of_weights

    # --- Output the results as requested ---
    print("The Crawley-Nordström invariant (c_N) is given by the formula:")
    print("c_N = d - Σwᵢ")
    print(f"\nThe given weights are: {weights}")
    
    # Format the sum of weights for the equation
    weights_str = " + ".join(map(str, weights))
    print(f"The sum of the weights is: Σwᵢ = {weights_str} = {sum_of_weights}")

    print(f"For a Calabi-Yau variety, the degree 'd' equals the sum of the weights, so d = {d}.")

    print("\nPlugging these values into the formula:")
    # The final equation showing each number
    print(f"c_N = {d} - {sum_of_weights} = {c_N}")

    # This part is for the final answer format, but we'll print the value
    # to be captured by the system.
    # The final output format <<<...>>> will be added at the end of the response.
    # To avoid printing it here, we can store it.
    # For this exercise, let's just print the number that will go inside.
    # print(c_N, file=sys.stderr) # Example of printing for a different stream

if __name__ == "__main__":
    solve_crawley_nordstrom()