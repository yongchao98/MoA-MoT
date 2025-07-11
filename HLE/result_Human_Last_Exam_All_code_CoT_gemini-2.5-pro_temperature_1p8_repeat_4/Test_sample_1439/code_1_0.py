def solve_exponent_order():
    """
    Determines and explains the order in the coupling constant 'u' where the
    critical exponent 'nu' first receives a non-vanishing contribution in phi^4 theory.
    """
    # In mean-field theory (equivalent to the non-interacting case where u=0),
    # the critical exponent nu takes the value 1/2.
    nu_mean_field = 0.5
    
    # In the perturbative framework, we expand nu as a power series in the coupling constant 'u'.
    # The general form of the expansion is:
    # nu(u) = nu_0 + c_1*u^1 + c_2*u^2 + ...
    # where nu_0 is the mean-field value.
    
    # One-loop calculations in the renormalization group analysis show that the
    # coefficient c_1 is non-zero. This first correction arises from the Feynman
    # diagram that introduces a single interaction vertex, making the contribution
    # linear in the coupling constant 'u'.
    
    # Therefore, the full expression for nu begins as:
    # nu(u) = 1/2 + C * u^1 + O(u^2)
    # where C is a non-zero constant.
    
    order_of_first_contribution = 1
    
    print("The critical exponent nu can be expressed as a series in the coupling constant u.")
    print("The mean-field value (at u=0) is nu_0 = 1/2.")
    print("The expansion takes the form:")
    
    # Printing the equation as requested
    print(f"nu(u) = {nu_mean_field} + C * u^{order_of_first_contribution} + ... (where C is a non-zero constant)")
    
    print("\nThe initial non-vanishing contribution beyond the mean-field value comes from the term linear in u.")
    print(f"Thus, the specific order is: {order_of_first_contribution}")

solve_exponent_order()