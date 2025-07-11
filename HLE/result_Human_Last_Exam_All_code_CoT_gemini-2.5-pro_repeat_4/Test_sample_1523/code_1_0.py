import sympy

def display_inner_product_expression():
    """
    This function uses the sympy library to symbolically represent and explain
    the inner product (ϕ, D_ϕ) for a neutral scalar field in finite-temperature
    field theory.
    """

    # --- Symbolic Definitions ---
    # Define symbols for coordinates, mass, and inverse temperature
    tau, m, beta = sympy.symbols('τ m β', real=True, positive=True)
    d = sympy.Symbol('d', integer=True, positive=True)

    # Define x as a placeholder for the d spatial coordinates
    # We cannot represent a vector directly in this simple setup
    x = sympy.Symbol('x') 
    
    # Define ϕ as a function of the spatial coordinates and imaginary time
    phi = sympy.Function('ϕ')(x, tau)

    # --- The Operator D ---
    # The operator D is the Euclidean Klein-Gordon operator: D = -∂_τ² - ∇² + m²
    # We represent the Laplacian ∇²ϕ symbolically as a single term for clarity.
    laplacian_phi = sympy.Symbol('(∇²ϕ)')
    
    # The term Dϕ is the operator D acting on the field ϕ
    D_phi = -sympy.Derivative(phi, (tau, 2)) - laplacian_phi + m**2 * phi

    # --- The Inner Product (ϕ, D_ϕ) ---
    # The inner product is the integral of ϕ * (Dϕ) over spacetime.
    # The integrand is ϕ multiplied by Dϕ.
    integrand = phi * D_phi
    
    # --- Outputting the Explanation ---
    print("In finite-temperature field theory, the Euclidean action S[ϕ] for a free neutral scalar field is:")
    print("S[ϕ] = ∫ d^d x ∫_0^β dτ [ (1/2)(∂_τ ϕ)² + (1/2)(∇ϕ)² + (1/2)m²ϕ² ]")
    print("\nAfter integration by parts, this can be written as S[ϕ] = (1/2) (ϕ, D_ϕ),")
    print("where D is the Euclidean Klein-Gordon operator and ( , ) denotes the spacetime integral inner product.")
    
    print("\nThe inner product (ϕ, D_ϕ) is therefore given by the integral of the following expression:")
    print("(ϕ, D_ϕ) = ∫ d^d x ∫_0^β dτ [ ϕ(x,τ) * (Dϕ)(x,τ) ]")

    print("\nLet's break down the integrand: ϕ * (Dϕ)")
    print("--------------------------------------------------")
    
    # Print each term and its coefficient as requested
    
    print("1. Term from the second time derivative (-∂_τ²):")
    term1 = phi * (-sympy.Derivative(phi, (tau, 2)))
    sympy.pprint(term1, use_unicode=True)
    print("The numerical coefficient here is -1.")
    
    print("\n2. Term from the Laplacian (-∇²):")
    term2 = phi * (-laplacian_phi)
    sympy.pprint(term2, use_unicode=True)
    print("The numerical coefficient here is -1.")

    print("\n3. Term from the mass (m²):")
    term3 = phi * (m**2 * phi)
    sympy.pprint(term3, use_unicode=True)
    print("The numerical coefficient here is 1.")

    print("\n--------------------------------------------------")
    print("The complete expression for the integrand of (ϕ, D_ϕ) is the sum of these three terms:")
    sympy.pprint(integrand, use_unicode=True)
    print("\nThis integrand, integrated over Euclidean spacetime, gives the inner product (ϕ, D_ϕ).")
    print("This value is also equal to twice the Euclidean action: (ϕ, D_ϕ) = 2S[ϕ].")

if __name__ == '__main__':
    display_inner_product_expression()