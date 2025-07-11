import math

# Step 1: Define the given values from the problem description.
fundamental_loss_mw = 2.0  # The fundamental power loss in MW.
harmonic_loss_percentage = 0.10  # 10% additional loss for each harmonic order.

# Step 2: Calculate the power loss for each harmonic.
# The loss for the 3rd harmonic is 10% of the fundamental loss.
loss_3rd_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage

# The loss for the 5th harmonic is also 10% of the fundamental loss.
loss_5th_harmonic_mw = fundamental_loss_mw * harmonic_loss_percentage

# Step 3: Calculate the total power loss.
# The total loss is the sum of the fundamental loss and the losses from the harmonics.
total_loss_mw = fundamental_loss_mw + loss_3rd_harmonic_mw + loss_5th_harmonic_mw

# Step 4: Print the detailed calculation and the final result.
print("Calculation of Total Resistive Power Loss:")
print(f"Fundamental Loss = {fundamental_loss_mw} MW")
print(f"3rd Harmonic Loss (10% of Fundamental) = {loss_3rd_harmonic_mw:.2f} MW")
print(f"5th Harmonic Loss (10% of Fundamental) = {loss_5th_harmonic_mw:.2f} MW")
print("\nFinal Equation for Total Power Loss:")
print(f"Total Loss = Fundamental Loss + 3rd Harmonic Loss + 5th Harmonic Loss")
print(f"Total Loss = {fundamental_loss_mw} MW + {loss_3rd_harmonic_mw:.2f} MW + {loss_5th_harmonic_mw:.2f} MW")
print(f"Total Loss = {total_loss_mw:.2f} MW")

# Output the final numerical answer in the specified format.
# We are asked to calculate the total resistive power losses including the harmonic effects.
final_answer = total_loss_mw
print(f"\n<<<>>>\n{final_answer}")