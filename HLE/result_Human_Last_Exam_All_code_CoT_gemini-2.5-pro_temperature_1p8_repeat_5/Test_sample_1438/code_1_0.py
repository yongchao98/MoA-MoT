def berezin_integral(coefficients):
    """
    Simulates the Berezin integral for a function of a single Grassmann variable.
    A function f(psi) is represented by its Taylor expansion f(psi) = a + b*psi.
    Since psi*psi = 0, there are no higher order terms.
    The coefficients are provided as a list [a, b].

    According to Berezin integration rules:
    ∫ d(psi) f(psi) = ∫ d(psi) (a + b*psi) = a*∫d(psi) + b*∫d(psi)*psi = a*0 + b*1 = b.
    The integral simply returns the coefficient of the linear term.
    """
    if not isinstance(coefficients, list) or len(coefficients) != 2:
        raise TypeError("Input must be a list of two coefficients [a, b].")
    
    # The result of the integral is the coefficient 'b'.
    a, b = coefficients
    return b

def demonstrate_rules():
    """
    Demonstrates the fundamental rules of Berezin integration using the simulator.
    """
    print("In fermionic path integrals, the Pauli exclusion principle is enforced using Grassmann variables (e.g., ψ), where ψ²=0.")
    print("The integral over such variables is defined by the following rules:")
    print("-" * 60)

    # Rule 1: Integral of a constant (f(ψ) = 1, so a=1, b=0)
    constant_func_coeffs = [1, 0]
    result_1 = berezin_integral(constant_func_coeffs)
    print("1. The integral of a constant function f(ψ) = 1 is:")
    print(f"   ∫ dψ * {constant_func_coeffs[0]} = {result_1}")
    print("-" * 60)

    # Rule 2: Integral of the variable itself (f(ψ) = ψ, so a=0, b=1)
    linear_func_coeffs = [0, 1]
    result_2 = berezin_integral(linear_func_coeffs)
    print("2. The integral of the Grassmann variable itself f(ψ) = ψ is:")
    print(f"   ∫ dψ * ψ = {result_2}")
    print("-" * 60)
    
    print("These two rules, particularly the normalization to 1, define the measure for the integral.")

if __name__ == "__main__":
    demonstrate_rules()
