def calculate_critical_exponent_nu():
    """
    Calculates the critical exponent nu for a G₄-theoretical framework.

    The "G₄-theoretical framework" refers to a φ⁴ theory at its upper
    critical dimension, d=4. In this case, mean-field theory is exact,
    and the critical exponents take their mean-field values.

    We can also see this from the epsilon expansion for nu, where ε = 4 - d:
    ν = 1/2 + (1/12)ε + O(ε²)

    For d=4, ε = 0, and the higher-order terms vanish.
    """
    
    # Spatial dimension, inferred from the "G₄" notation
    d = 4
    
    # The parameter for the epsilon expansion
    epsilon = 4 - d
    
    # Coefficients from the epsilon expansion for ν
    # ν = term_0 + term_1 * ε + ...
    term_0_num = 1
    term_0_den = 2
    term_1_num = 1
    term_1_den = 12
    
    # Calculate the value of nu
    nu = (term_0_num / term_0_den) + (term_1_num / term_1_den) * epsilon
    
    # Print the explanation and the final equation
    print(f"The G₄-theoretical framework implies a spatial dimension d = {d}.")
    print(f"The epsilon expansion parameter is ε = 4 - d = {4} - {d} = {epsilon}.")
    print("The epsilon expansion for the critical exponent ν is: ν = 1/2 + (1/12)ε + ...")
    print("\nSubstituting ε = 0 into the equation:")
    print(f"ν = {term_0_num}/{term_0_den} + ({term_1_num}/{term_1_den}) * {epsilon}")
    print(f"ν = {term_0_num/term_0_den} + {0.0}")
    print(f"ν = {nu}")

if __name__ == "__main__":
    calculate_critical_exponent_nu()
