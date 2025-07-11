import sympy

def display_partition_function():
    """
    This function derives and prints the partition function Z
    for the given Hamiltonian in the grand canonical ensemble.
    """
    # Define symbols for the variables in the equation
    Z = sympy.Symbol('Z')
    beta = sympy.Symbol('beta')
    mu = sympy.Symbol('mu')

    # The numbers in the equation as requested
    num_one = 1
    num_two = 2

    # The expression for the partition function Z for a single bosonic mode.
    # The derivation is as follows:
    # 1. Z = Tr[exp(-beta * (H - mu*N))]
    # 2. Given H = -mu*N, the expression becomes Z = Tr[exp(-beta * (-mu*N - mu*N))]
    # 3. This simplifies to Z = Tr[exp(2*beta*mu*N)]
    # 4. For a single bosonic mode, the trace is a sum over particle numbers n = 0, 1, 2, ...
    #    Z = Sum_{n=0 to inf} [exp(2*beta*mu)]^n
    # 5. This is a geometric series which sums to 1 / (1 - x), where x = exp(2*beta*mu)
    
    # Constructing the final equation symbolically
    # Note: We use sympy's exp for symbolic representation
    expression = num_one / (num_one - sympy.exp(num_two * beta * mu))
    
    equation = sympy.Eq(Z, expression)

    # Print the final result
    print("The partition function Z for a system with the specified Hamiltonian is:")
    # We use pretty print for a nicer mathematical output
    sympy.pprint(equation, use_unicode=True)
    
    print("\nWhere:")
    print(f"  Z: The partition function")
    print(f"  {num_one}: The number one")
    print(f"  {num_two}: The number two")
    print(f"  beta: Inverse temperature (1/kT)")
    print(f"  mu: Chemical potential")
    print(f"  exp: The exponential function")
    print("\nThis result corresponds to a single bosonic mode and requires mu < 0 for convergence.")

if __name__ == '__main__':
    display_partition_function()