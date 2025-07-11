# Step 1: Define the given values
fundamental_loss_mw = 2.0  # Fundamental power loss in MW
harmonic_loss_factor = 0.10  # Additional loss per harmonic order (10%)

# Step 2: Calculate the loss due to the 3rd harmonic
loss_3rd_harmonic = fundamental_loss_mw * harmonic_loss_factor

# Step 3: Calculate the loss due to the 5th harmonic
loss_5th_harmonic = fundamental_loss_mw * harmonic_loss_factor

# Step 4: Calculate the total power loss
total_loss_mw = fundamental_loss_mw + loss_3rd_harmonic + loss_5th_harmonic

# Step 5: Print the final equation and the result
print("Calculation of Total Resistive Power Loss:")
print(f"Fundamental Loss = {fundamental_loss_mw} MW")
print(f"Additional Loss from 3rd Harmonic (10% of fundamental) = {loss_3rd_harmonic:.2f} MW")
print(f"Additional Loss from 5th Harmonic (10% of fundamental) = {loss_5th_harmonic:.2f} MW")
print("\nFinal Equation:")
print(f"Total Loss = {fundamental_loss_mw} MW + {loss_3rd_harmonic:.2f} MW + {loss_5th_harmonic:.2f} MW = {total_loss_mw:.2f} MW")

# The final answer is the total calculated power loss.
# The format <<<value>>> is used to indicate the final answer as requested.
print("\n<<<2.4>>>")