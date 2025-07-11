import sys

def get_shapley_formula(n):
    """
    For a given number of people n, this function calculates the coefficients
    for the formula that gives the Shapley value c_k for person p_k.
    """
    
    # Calculate the sum of the first n integers, S
    S = n * (n + 1) // 2
    
    # Calculate the sum of the first n squares, Q
    Q = n * (n + 1) * (2 * n + 1) // 6
    
    # The formula for c_k is derived as: c_k = A*k - B*k^2
    # where A = S^3 + S*Q and B = S^2
    
    # Calculate coefficient A for the k term
    A = S**3 + S * Q
    
    # Calculate coefficient B for the k^2 term
    B = S**2
    
    print(f"For n = {n}, the formula for the amount c_k for person p_k is:")
    
    # The final equation requires printing each number
    print(f"c_k = ({S}^3 + {S}*{Q}) * k - ({S}^2) * k^2")
    print(f"c_k = ({S**3} + {S*Q}) * k - {S**2} * k^2")
    print(f"c_k = {A} * k - {B} * k^2")


if __name__ == '__main__':
    # You can provide n as a command-line argument.
    # Example: python your_script_name.py 10
    try:
        n_val = int(sys.argv[1])
        if n_val <= 1:
            raise ValueError("n must be greater than 1.")
    except (IndexError, ValueError) as e:
        print(f"Usage: python <script> n, where n is an integer > 1.")
        print("Using a default value of n=4 for demonstration.")
        n_val = 4
    
    get_shapley_formula(n_val)
