def print_kappa_definition():
    """
    Prints the definition of the variable kappa based on the derivation.
    """
    # The variable sigma_x_sq represents the variance of the input signals, sigma_x^2.
    # The other variables correspond to the parameters given in the model description.
    kappa_definition = "kappa = - ( (mu**2 + (rho * mu) / phi) * (tau_u + tau_v) ) / sigma_x_sq"
    print("The definition of kappa is:")
    print(kappa_definition)

print_kappa_definition()