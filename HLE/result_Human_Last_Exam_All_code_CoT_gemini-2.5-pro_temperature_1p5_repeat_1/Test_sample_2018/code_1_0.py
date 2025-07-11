import math

def get_original_computational_factor():
    """
    Calculates and prints the value of the computational factor used in the 
    prior simulation-only work for the enthalpy-porosity method.
    """
    
    # The Carman-Kozeny source term added to the momentum equation is of the form:
    # S = -C * ((1-f)^2 / (f^3 + epsilon)) * u
    # where C is the computational factor (or mushy zone constant).

    # The prior, simulation-only work by Voller and Prakash (1987)
    # established a value for C. The question asks for this original value.
    base_value = 1.6
    exponent = 6

    # Calculate the final value
    computational_factor = base_value * (10 ** exponent)

    print("The computational factor 'C' in the prior simulation-only work was based on the following equation:")
    print(f"Value = {base_value} * 10^{exponent}")
    print("\nResult:")
    print(f"The originally used value for the computational factor was {computational_factor:e}")

if __name__ == '__main__':
    get_original_computational_factor()