import sys

def calculate_harmonic_losses():
    """
    Calculates the total resistive power loss including fundamental, 
    third, and fifth harmonic losses.
    """
    # Given values
    fundamental_loss_mw = 2.0
    harmonic_loss_percentage = 0.10  # 10%

    # Calculate additional losses for each harmonic
    # Loss from 3rd harmonic is 10% of the fundamental loss
    loss_3rd_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage
    
    # Loss from 5th harmonic is 10% of the fundamental loss
    loss_5th_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage

    # Calculate total power loss
    total_loss_mw = fundamental_loss_mw + loss_3rd_harmonic_mw + loss_5th_harmonic_mw

    # Print the explanation and the final equation
    print("The total resistive power loss is the sum of the fundamental loss and the additional losses from the 3rd and 5th harmonics.")
    print("\nCalculation:")
    print(f"Fundamental Loss = {fundamental_loss_mw} MW")
    print(f"Additional Loss from 3rd Harmonic = {loss_3rd_harmonic_mw:.1f} MW")
    print(f"Additional Loss from 5th Harmonic = {loss_5th_harmonic_mw:.1f} MW")
    
    print("\nFinal Equation:")
    # Using sys.stdout.write to prevent adding a newline, which can be useful in some contexts.
    # Here, we build the string first for clarity.
    equation_str = f"{fundamental_loss_mw} MW + {loss_3rd_harmonic_mw:.1f} MW + {loss_5th_harmonic_mw:.1f} MW = {total_loss_mw:.1f} MW"
    print(equation_str)
    
    # Returning the final numerical value for the answer block
    return total_loss_mw

if __name__ == '__main__':
    total_loss = calculate_harmonic_losses()
    # The final answer is formatted as requested.
    # The print statement below is for display purposes in a terminal. The actual answer format is handled separately.
    # print(f"\n<<< {total_loss} >>>")
    sys.stdout.write(f"\n<<<{total_loss:.1f}>>>")
