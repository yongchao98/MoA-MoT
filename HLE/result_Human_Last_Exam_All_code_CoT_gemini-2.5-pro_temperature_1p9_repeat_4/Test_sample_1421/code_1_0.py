def display_fermionic_partition_function():
    """
    Constructs and prints the formula for the fermionic partition function
    in the path integral formalism.
    """
    
    # Define the numbers present in the formula
    number_zero = 0
    number_two = 2

    # Build the formula string piece by piece for clarity
    partition_function = "Z"
    equals = "="
    functional_integral = "∫ Dψ̄ Dψ"
    exponential = "exp"
    open_bracket = "["
    close_bracket = "]"
    action_sign = "-"
    time_integral = f"∫_({number_zero} to β) dτ"
    space_integral = "∫ d^d x"
    field_left = "ψ̄(x, τ)"
    operator_part1 = " (∂_τ"
    operator_part2 = f"- (ħ² / ({number_two}*m))∇²"
    operator_part3 = "- μ)"
    field_right = "ψ(x, τ)"
    
    # Combine the pieces into the final equation string
    equation = (
        f"{partition_function} {equals} {functional_integral} {exponential} {open_bracket}"
        f"{action_sign} {time_integral} {space_integral} {field_left}"
        f"{operator_part1} {operator_part2} {operator_part3} {field_right} {close_bracket}"
    )

    print("The formula for the fermionic partition function (Z) is:")
    print(equation)
    print("\nThis formula is subject to the anti-periodic boundary condition in imaginary time:")
    print("ψ(x, 0) = -ψ(x, β)")

# Execute the function to print the formula
display_fermionic_partition_function()