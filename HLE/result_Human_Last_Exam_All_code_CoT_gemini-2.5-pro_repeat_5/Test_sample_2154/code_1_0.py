import math

def solve_picard_fuchs_order():
    """
    Calculates the minimal order u_r(n) of the Picard-Fuchs differential equation
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n = 3 to 12.
    """
    print("The problem is to find the minimal order u_r(n) of the Picard-Fuchs differential equation")
    print("for the period of the Hamiltonian system with potential V(q) = 1/2 * (q^2 - q^n).")
    print("This order is known to be equal to the genus of the associated hyperelliptic curve,")
    print("which is given by the formula: u_r(n) = floor((n - 1) / 2).")
    print("-" * 20)

    # Define the range for n
    n_values = range(3, 13)

    # A list to hold the results
    u_r_values = []

    print("Calculating the values for u_r(n) for n = 3, 4, ..., 12:")
    
    # Loop through each n, calculate u_r(n), and print the details
    for n in n_values:
        # In Python, integer division // on positive numbers is equivalent to floor()
        order = (n - 1) // 2
        u_r_values.append(order)
        # The prompt asks to "output each number in the final equation"
        print(f"u_r({n}) = floor(({n} - 1) / 2) = {order}")

    # Print the final result as a set expression
    result_string = ", ".join(map(str, u_r_values))
    print("\nThe complete set is:")
    print(f"{{u_r(3), u_r(4), ..., u_r(12)}} = {{{result_string}}}")

solve_picard_fuchs_order()