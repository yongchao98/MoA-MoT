class Blockchain:
    """
    A simple class to simulate a blockchain and demonstrate the concepts of
    Safety and Liveness.
    """
    def __init__(self):
        # The existing, confirmed blocks. We assume this part is immutable.
        self.chain = ["Block 0: Genesis"]
        self.pending_transactions = []
        # A flag to simulate a network stall.
        self.is_stalled = False

    def check_safety(self):
        """
        Checks the Safety property. In our simulation, safety means the
        original chain remains untampered. Since we don't modify the
        existing `self.chain` list (only append to it), safety is preserved.
        """
        print("Checking Safety: The existing chain is considered immutable and correct.")
        print("Result: Safety is NOT broken.")
        return True

    def check_liveness(self):
        """
        Checks the Liveness property. Liveness means the system can still
        add new blocks. This is directly affected by our `is_stalled` flag.
        """
        if self.is_stalled:
            print("Checking Liveness: The system is stalled and cannot add new blocks.")
            print("Result: Liveness IS broken.")
            return False
        else:
            print("Checking Liveness: The system is active and can make progress.")
            print("Result: Liveness is NOT broken.")
            return True

# --- Main Scenario ---
print("--- Scenario: No new blocks have been added for a long time. ---\n")

# 1. Create a blockchain instance.
my_blockchain = Blockchain()

# 2. A transaction is submitted, waiting to be included.
my_blockchain.pending_transactions.append("Alice sends 5 BTC to Bob")
print(f"A transaction is pending: {my_blockchain.pending_transactions[0]}\n")

# 3. A network-wide stall occurs (e.g., consensus failure, no miners).
# This prevents new blocks from being created.
print("*** A network stall occurs! ***\n")
my_blockchain.is_stalled = True

# 4. Evaluate the properties of the blockchain in this stalled state.
print("--- Evaluating Blockchain Properties ---")
safety_ok = my_blockchain.check_safety()
print("-" * 20)
liveness_ok = my_blockchain.check_liveness()
print("-" * 20)

print("\n--- Conclusion ---")
print("The simulation shows that even when the system stops producing blocks (a liveness failure),")
print("the integrity of the existing chain is not necessarily compromised (no safety failure).")
print("Therefore, it is for sure that liveness is broken, but safety may not be broken.")