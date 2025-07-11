def display_transformation_on_x():
    """
    This function prints the general representation of the Lie group transformation
    for the x-coordinate of the given PDE.
    The transformation depends on three arbitrary constants C1, C2, C3 from the
    Lie algebra, and the group parameter epsilon.
    """
    
    # Header
    print("The equation for the heat equation with a logarithmic source term is:")
    print("u_t = u_{xx} + (k_1*ln(u) + k_2)*u\n")
    
    # Infinitesimals
    print("The infinitesimal generators for the transformation are:")
    print("tau = C1")
    print("xi = C2 + C3 * exp(k1*t)")
    print("eta = (C4 - (1/2)*k1*C3*x) * exp(k1*t) * u")
    print("where C1, C2, C3, C4 are arbitrary constants.\n")
    
    # Finite transformation introduction
    print("The finite transformation x' for the x-coordinate is found by integrating:")
    print("dx'/d(epsilon) = xi(t'), with dt'/d(epsilon) = tau(t')")
    print("The result depends on whether the constant C1 (from tau) is zero.\n")
    
    # Case 1: C1 is not zero
    print("Case 1: C1 is not equal to 0")
    print("The transformation is:")
    # Print the equation term by term as requested
    print("x' = ", end="")
    print("x", end="")
    print(" + ", end="")
    print("C2*epsilon", end="")
    print(" + ", end="")
    print("(C3 / (k1*C1)) * exp(k1*t) * (exp(k1*C1*epsilon) - 1)")
    print("-" * 50)
    
    # Case 2: C1 is zero
    print("Case 2: C1 is equal to 0")
    print("The transformation is:")
    # Print the equation term by term
    print("x' = ", end="")
    print("x", end="")
    print(" + ", end="")
    print("(C2 + C3 * exp(k1*t)) * epsilon")
    print("-" * 50)

if __name__ == '__main__':
    display_transformation_on_x()
