import hashlib
import datetime

class Block:
    """
    A class representing a single block in the blockchain.
    Each block stores data for one software update.
    """
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data  # e.g., {'version': '1.1', 'file_hash': 'abc...'}
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """
        Calculates the SHA-256 hash of the block's contents.
        """
        # The 'equation' to create the hash is the concatenation of its properties.
        # We output the index as a key 'number' in the process.
        print(f"Hashing data for Block #{self.index}...")
        block_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(block_string.encode()).hexdigest()

class Blockchain:
    """
    A class to manage the chain of blocks.
    """
    def __init__(self):
        # The first block in the chain is the "Genesis Block"
        self.chain = [self.create_genesis_block()]

    def create_genesis_block(self):
        """
        Manually creates the first block with index 0.
        """
        return Block(0, datetime.datetime.now(), "Genesis Block: DNS Software v1.0", "0")

    def get_latest_block(self):
        """
        Returns the most recent block in the chain.
        """
        return self.chain[-1]

    def add_block(self, new_block):
        """
        Adds a new block to the chain after verifying it.
        """
        new_block.previous_hash = self.get_latest_block().hash
        new_block.hash = new_block.calculate_hash()
        self.chain.append(new_block)

    def is_chain_valid(self):
        """
        Verifies the integrity of the entire blockchain.
        It checks if hashes and links between blocks are correct.
        """
        print("\n--- Verifying Blockchain Integrity ---")
        for i in range(1, len(self.chain)):
            current_block = self.chain[i]
            previous_block = self.chain[i-1]

            # Check if the stored hash of the block is actually correct
            if current_block.hash != current_block.calculate_hash():
                print(f"Verification FAILED: The hash of Block #{current_block.index} is invalid.")
                return False

            # Check if the block points to the correct previous block
            if current_block.previous_hash != previous_block.hash:
                print(f"Verification FAILED: The chain is broken at Block #{current_block.index}.")
                return False
        print("Verification SUCCESS: The blockchain is valid and has not been tampered with.")
        return True

# --- Simulation ---

# 1. Create the secure ledger for our DNS software updates.
dns_update_ledger = Blockchain()
print("Created a new secure ledger for DNS software updates.")

# 2. Legitimate developers release new, signed updates.
print("\n--- Releasing Legitimate Updates ---")
update_v1_1 = {'version': '1.1', 'file_hash': 'a1b2c3d4e5f6'}
dns_update_ledger.add_block(Block(1, datetime.datetime.now(), update_v1_1, ""))

update_v1_2 = {'version': '1.2', 'file_hash': 'f6e5d4c3b2a1'}
dns_update_ledger.add_block(Block(2, datetime.datetime.now(), update_v1_2, ""))

print("\n--- Current Valid Blockchain ---")
for block in dns_update_ledger.chain:
    print(f"Index: {block.index}")
    print(f"Data: {block.data}")
    print(f"Hash: {block.hash}")
    print(f"Previous Hash: {block.previous_hash}\n")

# 3. Verify the chain. It should be valid.
dns_update_ledger.is_chain_valid()

# 4. An adversary attempts to tamper with a previous update to inject a malicious implant.
print("\n--- Adversary Tampering Attempt ---")
print("Adversary modifies the data in Block #1 to point to a malicious file...")
# The hash in the data would change if the file changed.
malicious_data = {'version': '1.1', 'file_hash': 'MALICIOUS_HASH_HERE_ca978112ca1bbdc'}
dns_update_ledger.chain[1].data = malicious_data

# 5. Run the verification process again.
# The system now detects the tampering because the hashes no longer match.
# The 'equation' for calculating the hash now produces a different result for Block 1,
# which breaks the entire chain's integrity.
dns_update_ledger.is_chain_valid()

<<<G>>>