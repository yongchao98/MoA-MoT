def display_fermionic_partition_function():
    """
    This function constructs and prints the formula for the fermionic
    partition function in the imaginary time path integral formalism.
    """

    # --- Step 1: Define the mathematical symbols using Unicode ---
    integral = "\u222B"
    partial_deriv = "\u2202"
    psi = "\u03C8"
    psi_star = f"{psi}*"
    tau = "\u03C4"
    beta = "\u03B2"

    # --- Step 2: Define the components of the formula as strings ---
    
    # The partition function symbol
    partition_function_Z = "Z"
    
    # The functional integration measure over Grassmann fields
    measure = f"D{psi_star} D{psi}"
    
    # The kinetic and Hamiltonian parts of the Euclidean Lagrangian
    kinetic_term = f"{psi_star}({tau}) {partial_deriv}{tau} {psi}({tau})"
    hamiltonian_term = f"H({psi_star}({tau}), {psi}({tau}))"
    
    # The full Euclidean action, S, which is the integral of the Lagrangian
    action_integral = f"{integral}_0^{beta} d{tau} [ {kinetic_term} + {hamiltonian_term} ]"
    
    # --- Step 3: Assemble and print the full formula with explanations ---

    print("The formula for the fermionic partition function Z in imaginary time is given by a path integral over anticommuting Grassmann fields.")
    print("-" * 80)
    print("Explanation of terms:")
    print(f"  Z          : The partition function, Tr(exp(-{beta}H)).")
    print(f"  {integral} {measure} : The functional integral over all anti-periodic field configurations.")
    print(f"  {beta}          : Inverse temperature (1 / k_B T).")
    print(f"  {tau}          : Imaginary time.")
    print(f"  {psi}({tau}), {psi_star}({tau}) : Fermionic (Grassmann) fields.")
    print(f"  H({psi_star}, {psi})  : The classical Hamiltonian with operators replaced by Grassmann fields.")
    print(f"  S          : The Euclidean action, given by the integral: S = {action_integral}")
    print("\nCrucial Boundary Condition:")
    print(f"The fields must be anti-periodic in imaginary time: {psi}({beta}) = -{psi}(0)")
    print("-" * 80)
    
    print("\nThe final formula, showing each component, is:")

    # Using an f-string to build the final equation from its components
    # This fulfills the requirement to output each part of the final equation.
    final_equation = (
        f"{partition_function_Z} = "
        f"{integral} {measure} "
        f"exp( - ( {action_integral} ) )"
    )
    
    print(final_equation)


if __name__ == "__main__":
    display_fermionic_partition_function()