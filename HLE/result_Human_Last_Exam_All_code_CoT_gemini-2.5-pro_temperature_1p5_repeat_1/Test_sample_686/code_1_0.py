def solve_ising_susceptibility():
    """
    This function provides the derived expression for the magnetic susceptibility
    of an Ising model on a sparse random graph in the paramagnetic phase.
    It prints the final equation and its components.
    """
    # The final expression is derived as chi = (beta * c * tanh(beta*J)) / (1 - (c-1)*tanh(beta*J))
    
    chi = "chi"
    beta = "beta"
    c = "c"
    J = "J"
    tanh = "tanh"
    one = "1"
    
    numerator = f"{beta} * {c} * {tanh}({beta} * {J})"
    denominator = f"{one} - ({c} - {one}) * {tanh}({beta} * {J})"
    
    print("The final expression for the magnetic susceptibility is:")
    print(f"{chi} = ({numerator}) / ({denominator})")
    
    print("\nHere is a breakdown of each component in the final equation:")
    print(f"{chi}: The magnetic susceptibility.")
    print(f"{beta}: The inverse temperature (often 1/kT).")
    print(f"{c}: The connectivity (average degree) of the graph.")
    print(f"{J}: The homogeneous coupling constant.")
    print(f"{tanh}: The hyperbolic tangent function.")
    print(f"{one}: The number one.")

solve_ising_susceptibility()
<<<chi = (beta * c * tanh(beta * J)) / (1 - (c - 1) * tanh(beta * J))>>>