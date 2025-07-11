import sys

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau Link.
    """
    # The weights of the ambient space
    weights = [22, 29, 49, 50, 75]

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # For a Calabi-Yau manifold, the degree 'd' of the defining polynomial
    # is equal to the sum of the weights.
    d = sum_of_weights

    # The Crawley-Nordström invariant 'c' is defined as c = d - sum(w_i).
    invariant = d - sum_of_weights
    
    # Print the final equation as requested, showing each number.
    print(f"The weights are: {weights}")
    weights_sum_str = " + ".join(map(str, weights))
    print(f"The sum of the weights is: {weights_sum_str} = {sum_of_weights}")
    print(f"For a Calabi-Yau Link, the degree 'd' of the polynomial must equal the sum of the weights, so d = {d}.")
    print("\nThe Crawley-Nordström invariant 'c' is calculated as: c = d - (sum of weights)")
    print(f"c = {d} - {sum_of_weights}")
    print(f"c = {invariant}")

solve_crawley_nordstrom()
<<<0>>>