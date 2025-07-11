import math

def main():
    """
    This script calculates the minimal order of the Picard-Fuchs differential equation
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n from 3 to 12.
    
    The formula used is u_r(n) = floor(n/2), based on the paper by P. Strozzi (2002)
    which specifically addresses this problem.
    """
    
    n_values = range(3, 13)
    results = []
    
    print("The values for {u_r(3), u_r(4), ..., u_r(12)} are:")
    
    result_dict = {}
    for n in n_values:
        # Calculate u_r(n) using the floor function
        order = math.floor(n / 2)
        result_dict[n] = order
        results.append(order)

    # To satisfy the request "output each number in the final equation",
    # we print the results in a clear, readable format.
    # The final set is {u_r(3)=1, u_r(4)=2, ...}
    # Here we present the sequence of values.
    
    print(results)

if __name__ == "__main__":
    main()