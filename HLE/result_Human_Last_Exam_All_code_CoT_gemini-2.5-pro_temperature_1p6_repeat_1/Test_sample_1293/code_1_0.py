import hashlib
import datetime as date

class Block:
    """A block in the blockchain, representing a software update record."""
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data  # e.g., 'DNS Server v1.1, file_hash: 5f8d...'
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """Calculates the SHA-256 hash of the block's contents."""
        hash_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(hash_string.encode()).hexdigest()

class Blockchain:
    """A chain of blocks to securely track software updates."""
    def __init__(self):
        self.chain = [self.create_genesis_block()]

    def create_genesis_block(self):
        """Creates the very first block in the chain."""
        return Block(0, date.datetime.now(), "Genesis Block: Initial DNS Software v1.0", "0")

    def get_latest_block(self):
        """Returns the most recent block in the chain."""
        return self.chain[-1]

    def add_block(self, new_block):
        """Adds a new block (update) to the chain after verification."""
        new_block.previous_hash = self.get_latest_block().hash
        new_block.hash = new_block.calculate_hash()
        self.chain.append(new_block)

    def is_chain_valid(self):
        """Verifies the integrity of the entire blockchain."""
        for i in range(1, len(self.chain)):
            current_block = self.chain[i]
            previous_block = self.chain[i-1]

            # Check 1: Check if the hash of the block is correct
            if current_block.hash != current_block.calculate_hash():
                print(f"\nVerification FAILED: Current hash does not match calculated hash for Block {current_block.index}.")
                return False

            # Check 2: Check if the block points to the correct previous block
            if current_block.previous_hash != previous_block.hash:
                print(f"\nVerification FAILED: Chain broken. Previous hash mismatch at Block {current_block.index}.")
                return False
        
        print("\nVerification SUCCESS: Blockchain integrity is intact.")
        return True

# --- Simulation ---

# 1. Initialize the secure update ledger
update_ledger = Blockchain()
print("A secure blockchain ledger for DNS software updates has been created.")

# 2. Release a new legitimate software update
print("Releasing DNS Server v1.1...")
update_v1_1 = Block(1, date.datetime.now(), {"version": "1.1", "notes": "Security patch for CVE-2023-1234"}, "")
update_ledger.add_block(update_v1_1)
print("Update v1.1 successfully added to the ledger.")

# 3. Release another legitimate update
print("Releasing DNS Server v1.2...")
update_v1_2 = Block(2, date.datetime.now(), {"version": "1.2", "notes": "Performance improvements"}, "")
update_ledger.add_block(update_v1_2)
print("Update v1.2 successfully added to the ledger.")

# 4. Verify the chain - it should be valid
print("\n--- Running verification on the legitimate update chain... ---")
update_ledger.is_chain_valid()
print("Current chain state:")
for block in update_ledger.chain:
    print(f"  Index: {block.index}, Data: {block.data}, Hash: {block.hash[:10]}..., Prev Hash: {block.previous_hash[:10]}...")

# 5. Simulate an ATTACK: an adversary tampers with a previous update
print("\n--- !!! Adversary attempting a supply-chain attack !!! ---")
print("Tampering with Block 1 data to inject malicious code...")
# The attacker changes the data, which should invalidate all subsequent hashes
update_ledger.chain[1].data = {"version": "1.1-MALICIOUS", "notes": "Trojan implanted"}

# 6. Run verification again. The system should now detect the tampering.
print("\n--- Running verification after the tampering attempt... ---")
update_ledger.is_chain_valid()
