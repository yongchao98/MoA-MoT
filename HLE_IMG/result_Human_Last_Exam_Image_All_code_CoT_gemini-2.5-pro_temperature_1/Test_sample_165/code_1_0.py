def calculate_harmonic_losses():
    """
    Calculates the resistive power losses due to third and fifth harmonic distortions.
    """
    # 1. Define the given values
    fundamental_loss_mw = 2.0
    harmonic_loss_percentage = 0.10  # 10%

    # 2. Calculate the power loss for each harmonic
    # The loss for each harmonic is an additional 10% of the fundamental loss.
    loss_3rd_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage
    loss_5th_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage

    # 3. Calculate the total power loss due to harmonics
    total_harmonic_loss_mw = loss_3rd_harmonic_mw + loss_5th_harmonic_mw

    # 4. Print the calculation steps and the final result
    print("Problem: Calculate the resistive power losses due to 3rd and 5th harmonics.")
    print(f"Given Fundamental Power Loss = {fundamental_loss_mw} MW")
    print(f"Given Additional Loss per Harmonic = {harmonic_loss_percentage*100}% of Fundamental Loss\n")

    print("Step 1: Calculate loss due to the 3rd harmonic.")
    print(f"Loss_3rd = {fundamental_loss_mw} MW * {harmonic_loss_percentage} = {loss_3rd_harmonic_mw:.1f} MW\n")

    print("Step 2: Calculate loss due to the 5th harmonic.")
    print(f"Loss_5th = {fundamental_loss_mw} MW * {harmonic_loss_percentage} = {loss_5th_harmonic_mw:.1f} MW\n")

    print("Step 3: Calculate the total power loss from harmonics.")
    print("Total Harmonic Loss = Loss_3rd + Loss_5th")
    print(f"Total Harmonic Loss = {loss_3rd_harmonic_mw:.1f} MW + {loss_5th_harmonic_mw:.1f} MW = {total_harmonic_loss_mw:.1f} MW")

    # The final answer to be extracted
    # <<< {total_harmonic_loss_mw} >>>

if __name__ == "__main__":
    calculate_harmonic_losses()