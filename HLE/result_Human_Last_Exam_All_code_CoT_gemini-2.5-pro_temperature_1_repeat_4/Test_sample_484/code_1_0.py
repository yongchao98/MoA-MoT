def display_final_equation():
    """
    This function prints the derived mathematical expression for the synaptic efficacy dynamics
    and the definitions of the variables involved.
    """

    # Using unicode characters for Greek letters for better readability in terminals that support it.
    tau_w = "\u03c4_w"
    dw_dt = "dw\u1d62/dt"
    beta = "\u03b2"
    u_i = "u\u1d62"
    v_i = "v\u1d62"
    rho = "\u03c1"
    eta = "\u03b7"
    alpha = "\u03b1"
    phi = "\u03d5"

    # Constructing the final equation string
    lhs = f"{tau_w} {dw_dt}"
    rhs_numerator = f"{beta}{u_i}({v_i} - ({rho}(1 - {eta}) - {eta}))"
    rhs_denominator = f"1 + {v_i}"
    equation = f"{lhs} = ({rhs_numerator}) / ({rhs_denominator})"

    # Print the final result and definitions
    print("The derived expression for the change in synaptic efficacy is:")
    print(equation)
    print("\nEach term in the final equation is presented above.")
    
    print("\nThe variables in the reduced system are defined as:")
    print(f"w\u1d62: Synaptic efficacy.")
    print(f"{u_i} = \u03a3\u2c7c w\u2c7c x\u2c7c: Postsynaptic accumulator, representing the total weighted input (proportional to postsynaptic calcium \u0232).")
    print(f"{v_i} = {phi}x\u1d62: Presynaptic accumulator, representing presynaptic activity (proportional to presynaptic MMP9 M\u1d62).")
    print(f"{rho} = -{alpha}/{beta}: A constant representing the relative strength of LTD versus LTP.")
    print(f"The parameters {beta}, {eta}, and {phi} are constants from the original model.")

if __name__ == '__main__':
    display_final_equation()