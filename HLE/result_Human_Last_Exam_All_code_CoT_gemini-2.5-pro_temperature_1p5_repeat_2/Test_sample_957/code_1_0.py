import math

def print_final_equations():
    """
    This function prints the final equations for the reflection and transmission
    coefficients of a Fabry-Pérot cavity as given in the correct answer choice.
    The instruction 'output each number in the final equation' is interpreted
    as highlighting the numerical constants present in the formulas.
    """
    print("The derived transmission coefficient (τ) is:")
    
    # Equation for tau
    tau_numerator = "τ_m^2 * e^(i*k_0*d)"
    tau_denominator = "1 - ρ_m^2 * e^(i*2*k_0*d)"
    print(f"τ = ({tau_numerator}) / ({tau_denominator})\n")
    
    print("Numbers present in the transmission equation:")
    print("In the term τ_m^2, the number is: 2")
    print("In the term k_0, the number is: 0")
    print(f"In the term e^(i*k_0*d), the implied coefficient of d is: 1")
    print("In the denominator, the constant term is: 1")
    print("In the term ρ_m^2, the number is: 2")
    print("In the term e^(i*2*k_0*d), the coefficient of k_0*d is: 2")
    print("In the term k_0, the number is: 0")

    print("\n" + "="*50 + "\n")

    print("The corresponding reflection coefficient (ρ) from the selected option is:")

    # Equation for rho
    rho_numerator = "1 - (ρ_m - τ_m^2) * e^(i*2*k_0*d) * ρ_m"
    rho_denominator = "1 - ρ_m^2 * e^(i*2*k_0*d)"
    print(f"ρ = ({rho_numerator}) / ({rho_denominator})\n")

    print("Numbers present in the reflection equation:")
    print("In the numerator, the constant term is: 1")
    print("In the term τ_m^2, the number is: 2")
    print("In the term e^(i*2*k_0*d), the coefficient of k_0*d is: 2")
    print("In the term k_0, the number is: 0")
    print("In the denominator, the constant term is: 1")
    print("In the term ρ_m^2, the number is: 2")
    print("In the term e^(i*2*k_0*d), the coefficient of k_0*d is: 2")
    print("In the term k_0, the number is: 0")

if __name__ == '__main__':
    print_final_equations()