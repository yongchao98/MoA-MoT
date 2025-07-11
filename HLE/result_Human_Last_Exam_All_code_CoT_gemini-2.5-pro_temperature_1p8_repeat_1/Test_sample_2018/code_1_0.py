import math

def find_computational_factor():
    """
    This function identifies and presents the computational factor from the
    original enthalpy-porosity method paper by Voller and Prakash (1987).
    """

    # The value is given in the paper as a constant C.
    # The question asks for the value from the "prior published simulation-only work".
    # This refers to the 1987 paper by V.R. Voller and C. Prakash.
    
    # In Section 3.2 of Voller and Prakash (1987), the authors state:
    # "In the computations presented here C was set to 1.6 X 10^6..."
    
    base = 1.6
    exponent = 6
    
    # Calculate the full value of the constant
    computational_factor = base * (10**exponent)

    print("The method described is the enthalpy-porosity technique for modeling phase change.")
    print("The 'computational factor' is the mushy zone constant, C, in the momentum source term.")
    print("\nThe original value for this constant comes from the foundational simulation-only paper:")
    print("V. R. Voller and C. Prakash, Int. J. Heat Mass Transfer, Vol. 30, No. 8, pp. 1709-1719, 1987.")
    print("\nIn this paper, the authors specified the value used for their simulations.")
    print("\nThe equation for the computational factor is formed from the following numbers:")
    print(f"Base: {base}")
    print(f"Exponent: {exponent}")
    
    print("\nFinal Equation:")
    # The prompt requires outputting each number in the final equation.
    # We will format it as base * 10^exponent = value
    print(f"{base} * 10^{exponent} = {computational_factor:,.0f}")

find_computational_factor()
