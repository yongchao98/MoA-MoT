import cmath
import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the transmission lines
    based on the provided power system configuration and load conditions.
    """
    # Step 1: Define Constants and Parameters
    V_LL = 230e3  # Line-to-line voltage in Volts
    PF = 0.9      # Lagging power factor

    # Active power loads in Watts
    P_L_base = 50e6 / 3  # Base load at buses 5, 6, 8
    P_L_add_8 = 100e6    # Additional load at bus 8

    # Resistances of transmission lines in Ohms from the table
    R = {
        "T4_5": 5.29,
        "T4_6": 8.99,
        "T7_8": 4.50,
        "T8_9": 6.30,
        "T5_7": 16.93,
        "T6_9": 20.63
    }

    # Step 2: Calculate Complex Apparent Power for Loads
    # For a lagging PF, S = P + jQ. Q = P * tan(acos(PF)).
    tan_phi = math.tan(math.acos(PF))

    S_L5 = P_L_base * (1 + 1j * tan_phi)
    S_L6 = P_L_base * (1 + 1j * tan_phi)
    
    P_L8_total = P_L_base + P_L_add_8
    S_L8 = P_L8_total * (1 + 1j * tan_phi)

    # Step 3: Apply "Direct Feed" assumption for line flows
    # This model assumes loads are fed from the most direct generator path
    # and the load at Bus 8 is split between the two feeding lines.
    S_flow = {
        "T4_5": S_L5,
        "T4_6": S_L6,
        "T7_8": S_L8 / 2,
        "T8_9": S_L8 / 2,
        "T5_7": 0, # Assuming no flow in tie-lines for this simplified model
        "T6_9": 0  # Assuming no flow in tie-lines for this simplified model
    }

    # Step 4 & 5: Calculate power loss for each line
    P_loss = {}
    print("Calculating power loss for each transmission line (in MW):")

    for line, S in S_flow.items():
        S_mag = abs(S)
        R_line = R[line]
        # Using the formula for 3-phase power loss: P_loss = |S|^2 * R / V_LL^2
        loss_watts = (S_mag**2 * R_line) / (V_LL**2) if V_LL != 0 else 0
        P_loss[line] = loss_watts / 1e6 # Convert to MW

    # Step 6: Sum the losses
    total_loss_mw = sum(P_loss.values())

    # Step 7: Print the final calculation and result
    print("\nIndividual Line Losses (MW):")
    for line, loss in P_loss.items():
        print(f"  P_loss({line}) = {loss:.5f} MW")

    print("\nFinal Calculation of Total Loss:")
    loss_components = [f"{loss:.5f}" for loss in P_loss.values() if loss > 0]
    print("Total Loss = " + " + ".join(loss_components) + " MW")
    
    print(f"\nTotal power losses = {total_loss_mw:.3f} MW")


calculate_power_losses()
<<<0.950>>>