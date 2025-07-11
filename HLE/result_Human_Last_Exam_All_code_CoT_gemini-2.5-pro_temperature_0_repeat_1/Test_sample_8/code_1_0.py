def solve_conormal_space():
    """
    Calculates the resulting conormal space for R(sigma)f.
    """
    # The initial function f is in the conormal space A^{s}(X)
    # where s = 2 + alpha. We represent the constant part of the order.
    initial_order_k = 2
    alpha_symbol = "alpha"

    # The resolvent R(sigma) = (Box_g - sigma^2)^{-1} is the inverse of a
    # second-order operator. In the relevant pseudodifferential calculus,
    # this makes it an operator of order -2.
    operator_order = -2

    # An operator of order m maps A^s(X) to A^{s-m}(X).
    # So, the change in regularity order is -m.
    regularity_change = -operator_order

    # Calculate the final regularity order s' = s - m.
    final_order_k = initial_order_k - operator_order

    # Print the step-by-step derivation of the final space.
    print("Let the initial conormal regularity be s = 2 + alpha.")
    print("The resolvent R(sigma) is a pseudodifferential operator of order -2.")
    print("Applying an operator of order m = -2 maps A^s(X) to A^{s-m}(X).")
    print("The new regularity s' is calculated as follows:")
    
    # Output the final equation with each number
    print(f"s' = (s) - (m)")
    print(f"s' = ({initial_order_k} + {alpha_symbol}) - ({operator_order})")
    print(f"s' = {initial_order_k} + {alpha_symbol} + {regularity_change}")
    print(f"s' = {final_order_k} + {alpha_symbol}")

    print(f"\nTherefore, R(sigma)f belongs to the conormal space A^({final_order_k}+{alpha_symbol})(X).")

solve_conormal_space()