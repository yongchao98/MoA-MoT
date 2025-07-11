import sympy

def solve_minimal_surface_index():
    """
    Calculates the Morse index of a minimal surface given its Gauss map.

    The formula for the Morse index of a properly embedded minimal surface of
    genus zero is: Index = d_N + d_S - 1
    where:
    - d_N is the total number of poles of the Gauss map.
    - d_S is the total number of zeros of the Gauss map.
    """
    z = sympy.symbols('z')
    
    # Gauss map g(z) = P(z) / Q(z)
    P_str = "z"
    Q_str = "z**3 + 2"
    
    P = sympy.poly(P_str, z)
    Q = sympy.poly(Q_str, z)
    
    # d_N is the number of poles, which is the degree of the denominator Q(z)
    d_N = sympy.degree(Q)
    
    # d_S is the number of zeros, which must equal the degree of the map.
    # The degree of the map g(z) = P(z)/Q(z) is max(degree(P), degree(Q)).
    # Alternatively, we calculate it from finite zeros and behavior at infinity.
    
    # Number of finite zeros is the degree of the numerator P(z)
    num_finite_zeros = sympy.degree(P)
    
    # Order of the zero/pole at infinity is deg(Q) - deg(P)
    # If > 0, it's a zero at infinity. If < 0, it's a pole.
    order_at_infinity = sympy.degree(Q) - sympy.degree(P)
    
    if order_at_infinity < 0:
      # This case doesn't apply here but is included for completeness
      num_zeros_at_infinity = 0
    else:
      num_zeros_at_infinity = order_at_infinity

    # d_S is the total number of zeros
    d_S = num_finite_zeros + num_zeros_at_infinity
    
    # Morse Index calculation
    morse_index = d_N + d_S - 1
    
    print("For the Gauss map g(z) = z / (z^3 + 2):")
    print("-" * 40)
    print(f"The number of poles (d_N) is the degree of the denominator ({Q_str}), which is {d_N}.")
    print(f"The number of finite zeros is the degree of the numerator ({P_str}), which is {num_finite_zeros}.")
    print(f"The order of the zero at infinity is {num_zeros_at_infinity}.")
    print(f"The total number of zeros (d_S) is the sum of finite and infinite zeros: {num_finite_zeros} + {num_zeros_at_infinity} = {d_S}.")
    print("\nUsing the formula: Index = d_N + d_S - 1")
    print(f"Index = {d_N} + {d_S} - 1 = {morse_index}")
    print("-" * 40)
    print(f"The Morse index of the surface M is {morse_index}.")

solve_minimal_surface_index()