import sys

def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under selfish mining.

    Args:
        beta (float): The adversary's portion of mining power. Must be between 0 and 0.5.
        p (float): The probability that honest miners choose the adversary's block in a tie. Must be between 0 and 1.
    """
    # --- Parameter Validation ---
    if not 0 < beta < 0.5:
        print(f"Error: The adversary's mining power beta (β={beta}) must be between 0 and 0.5 for this model to be valid.", file=sys.stderr)
        print("If beta >= 0.5, the selfish miner's chain is expected to grow indefinitely, and the chain quality approaches 0.", file=sys.stderr)
        return
    if not 0 <= p <= 1:
        print(f"Error: The probability p (p={p}) must be between 0 and 1.", file=sys.stderr)
        return

    print(f"Calculating Chain Quality for β = {beta} and p = {p}\n")

    # --- Step 1: Calculate expected blocks in the "Honest Phase" ---
    # The honest phase lasts until the selfish miner finds a block.
    # The number of blocks in this phase follows a geometric distribution with success probability beta.
    # Expected number of trials is 1/beta. This includes the final selfish block.
    # So, the number of honest blocks preceding it is (1/beta - 1).
    exp_blocks_honest_phase = (1 / beta) - 1
    # In this phase, all blocks added to the chain are honest.
    exp_honest_blocks_honest_phase = exp_blocks_honest_phase
    
    # --- Step 2: Calculate expected blocks in the "Selfish Phase" ---
    # These formulas are derived from solving a system of linear equations
    # representing the state machine of the selfish mining attack.
    
    # E[H_s]: Expected number of HONEST blocks added during the selfish phase.
    exp_honest_blocks_selfish_phase = ((1 - beta)**2) * (2 - p)
    
    # E[L_s]: Expected TOTAL number of blocks added during the selfish phase.
    exp_total_blocks_selfish_phase = 2 + (beta**2) / (1 - beta)

    # --- Step 3: Calculate total expected blocks for a full cycle ---
    total_honest_blocks = exp_honest_blocks_honest_phase + exp_honest_blocks_selfish_phase
    total_chain_length = exp_blocks_honest_phase + exp_total_blocks_selfish_phase
    
    # --- Step 4: Calculate the final chain quality ---
    if total_chain_length == 0:
        chain_quality = 0
    else:
        chain_quality = total_honest_blocks / total_chain_length
        
    # --- Step 5: Print the results ---
    print("The final equation is in the form: (E[H_h] + E[H_s]) / (E[L_h] + E[L_s])\n")
    print(f"Expected Honest Blocks in Honest Phase (E[H_h]): {exp_honest_blocks_honest_phase}")
    print(f"Expected Honest Blocks in Selfish Phase (E[H_s]): {exp_honest_blocks_selfish_phase}")
    print(f"Expected Total Blocks in Honest Phase (E[L_h]): {exp_blocks_honest_phase}")
    print(f"Expected Total Blocks in Selfish Phase (E[L_s]): {exp_total_blocks_selfish_phase}")
    print("-" * 30)
    print("Final Equation with numbers:")
    print(f"({exp_honest_blocks_honest_phase:.4f} + {exp_honest_blocks_selfish_phase:.4f}) / ({exp_blocks_honest_phase:.4f} + {exp_total_blocks_selfish_phase:.4f})")
    print(f"= {total_honest_blocks:.4f} / {total_chain_length:.4f}")
    print(f"= {chain_quality:.4f}\n")
    
    # The problem asks for the answer to be returned in a specific format
    # at the very end. Let's return the final value in that format.
    print(f"The expected chain quality is: {chain_quality}")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    # Default values for demonstration if not provided via command line
    beta_val = 0.3
    p_val = 0.5
    
    if len(sys.argv) == 3:
        try:
            beta_val = float(sys.argv[1])
            p_val = float(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide numeric values for beta and p.", file=sys.stderr)
            sys.exit(1)
            
    calculate_chain_quality(beta_val, p_val)
