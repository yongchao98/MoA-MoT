def scalar_field_inner_product():
    """
    This function constructs and prints the equation for the inner product (ϕ, Dϕ)
    for a neutral scalar field in finite-temperature field theory.
    """
    # Define unicode symbols for mathematical notation
    integral = "\u222B"
    dtau = "d\u03c4"  # dτ
    ddx = "d\u1d48x"   # dᵈx (d-dimensional spatial volume element)
    partial = "\u2202"
    tau_sub = "\u208c" # Subscript τ
    nabla = "\u2207"
    phi = "\u03d5"
    sq = "\u00b2"      # Superscript 2

    # As requested, define the numbers (coefficients) to be included in the final equation.
    # For the standard action, these coefficients are all 1.
    coeff_kinetic_time = 1
    coeff_kinetic_space = 1
    coeff_mass = 1
    
    # Construct the individual terms of the integrand
    term1 = f"{coeff_kinetic_time}({partial}{tau_sub}{phi}){sq}"
    term2 = f"{coeff_kinetic_space}({nabla}{phi}){sq}"
    term3 = f"{coeff_mass}m{sq}{phi}{sq}"
    
    # Assemble the full equation as a string
    integrand = f"{term1} + {term2} + {term3}"
    full_equation = f"({phi}, D{phi}) = {integral} {dtau} {ddx} [ {integrand} ]"
    
    # Print the final result
    print(full_equation)

scalar_field_inner_product()