import hashlib
import datetime

class Block:
    """Represents a block in the blockchain, e.g., a software update record."""
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data # e.g., {'version': '1.2.3', 'file_hash': 'abc...', 'signature': '...'}
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """Calculates the SHA-256 hash of the block's contents."""
        block_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(block_string.encode()).hexdigest()

class Blockchain:
    """Manages the chain of blocks for software updates."""
    def __init__(self):
        self.chain = [self.create_genesis_block()]
        print("Blockchain for Secure Software Updates Initialized.")
        self.print_block(self.chain[0])

    def create_genesis_block(self):
        """Creates the very first block in the chain."""
        return Block(0, datetime.datetime.now(), "Genesis Block: Initial Software Release v1.0", "0")

    def get_latest_block(self):
        """Returns the most recent block in the chain."""
        return self.chain[-1]

    def add_block(self, new_block_data):
        """Adds a new valid block (update) to the chain."""
        latest_block = self.get_latest_block()
        new_block = Block(
            index=latest_block.index + 1,
            timestamp=datetime.datetime.now(),
            data=new_block_data,
            previous_hash=latest_block.hash
        )
        self.chain.append(new_block)
        print(f"\nNew Block #{new_block.index} Added to the Chain:")
        self.print_block(new_block)

    def print_block(self, block):
        """Prints the details of a single block."""
        print(f"  Index: {block.index}")
        print(f"  Timestamp: {block.timestamp}")
        print(f"  Data: {block.data}")
        print(f"  Previous Hash: {block.previous_hash}")
        print(f"  Hash: {block.hash}")

    def is_chain_valid(self):
        """Validates the integrity of the entire blockchain."""
        print("\n--- Running Integrity Check on the Blockchain ---")
        for i in range(1, len(self.chain)):
            current_block = self.chain[i]
            previous_block = self.chain[i-1]

            # Check if the stored hash of the block is correct
            if current_block.hash != current_block.calculate_hash():
                print(f"Validation Failed: Block #{current_block.index} has been tampered with. Calculated hash does not match stored hash.")
                return False
            
            # Check if the block points to the correct previous block
            if current_block.previous_hash != previous_block.hash:
                print(f"Validation Failed: Chain is broken at Block #{current_block.index}. Previous hash mismatch.")
                return False
        
        print("Integrity Check Passed: The blockchain is valid and secure.")
        return True

# --- Simulation ---
# 1. Create a secure update chain for a DNS software.
secure_update_chain = Blockchain()

# 2. Push a legitimate security patch.
patch_v1_1 = {
    'software': 'DNS Resolver',
    'version': '1.1',
    'description': 'Security patch for vulnerability CVE-2023-1001'
}
secure_update_chain.add_block(patch_v1_1)

# 3. Push a legitimate feature update.
update_v1_2 = {
    'software': 'DNS Resolver',
    'version': '1.2',
    'description': 'Added support for new record types'
}
secure_update_chain.add_block(update_v1_2)

# 4. Check the integrity of the chain.
secure_update_chain.is_chain_valid()

# 5. Simulate an adversary tampering with an old update to inject a malicious implant.
print("\n!!! Adversary attempts to tamper with Block #1 to inject malicious code. !!!")
secure_update_chain.chain[1].data = {
    'software': 'DNS Resolver',
    'version': '1.1-malicious',
    'description': 'Security patch for vulnerability CVE-2023-1001 with hidden backdoor'
}

# 6. Run the integrity check again. The system will now detect the tampering.
# The hash of block 1 will change, causing a mismatch with the 'previous_hash' stored in block 2.
secure_update_chain.is_chain_valid()
<<<G>>>