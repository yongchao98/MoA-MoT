import time

class SimpleBlockchain:
    """
    A simplified blockchain to demonstrate Safety vs. Liveness.
    """
    def __init__(self):
        # The chain of confirmed blocks (represents the safe, historical record)
        self.chain = ["Block 0 (Genesis)"]
        # Transactions waiting to be included in a block
        self.pending_transactions = []
        # A flag to simulate if the system is able to produce new blocks
        self.can_mine_blocks = True

    def check_liveness(self):
        """
        Liveness check: Are new transactions being processed?
        This is broken if the system can't mine new blocks while transactions are pending.
        """
        print("--- Checking Liveness ---")
        if self.pending_transactions and not self.can_mine_blocks:
            print("Result: LIVENESS IS BROKEN. Pending transactions exist, but no new blocks are being created.")
        else:
            print("Result: Liveness is OK. The system can process transactions.")

    def check_safety(self):
        """
        Safety check: Has the confirmed history of the chain been altered?
        In this simulation, we assume the chain's history is safe unless explicitly compromised.
        The scenario does not imply a compromise of past blocks.
        """
        print("--- Checking Safety ---")
        # In a real system, this would involve checking hashes and signatures.
        # Here, we just state that the history remains unaltered.
        print("Result: SAFETY IS NOT BROKEN. The integrity of the existing chain (Blocks 0, 1, etc.) is intact.")

# --- Simulation Setup ---
blockchain = SimpleBlockchain()
print("Blockchain created. Initially, the system is live.")
blockchain.pending_transactions.extend(["tx_A", "tx_B"])
blockchain.chain.append("Block 1, containing tx_A, tx_B")
blockchain.pending_transactions = []
print(f"Current Chain: {blockchain.chain}\n")


# --- The Problem Scenario Begins ---
print(">>> SCENARIO: System stalls. No new transactions have been included for 1 day. <<<\n")
# Simulate the system stall
blockchain.can_mine_blocks = False
# Users continue to submit new transactions
blockchain.pending_transactions.extend(["tx_C", "tx_D"])
print(f"State: System is stalled. can_mine_blocks = {blockchain.can_mine_blocks}")
print(f"State: New transactions are waiting. Pending Transactions: {blockchain.pending_transactions}\n")

# --- Analysis of the Scenario ---
# 1. Evaluate Liveness: The system has pending transactions but is unable to process them.
blockchain.check_liveness()

# 2. Evaluate Safety: The previously confirmed Block 0 and Block 1 are still secure and unchanged.
blockchain.check_safety()

print("\n--- Conclusion ---")
print("The analysis shows that it is for sure that liveness is broken, but safety may not be broken.")
