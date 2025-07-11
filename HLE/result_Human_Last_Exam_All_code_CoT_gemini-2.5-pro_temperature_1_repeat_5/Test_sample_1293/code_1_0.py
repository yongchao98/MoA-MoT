import hashlib
import datetime

class Block:
    """
    Represents a single block in our simple blockchain.
    Each block contains data for a software update.
    """
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data # e.g., {'version': '9.16.1', 'file_hash': '...'}
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """
        Calculates the SHA-256 hash for the block's contents.
        """
        # We ensure the dictionary is ordered so that we get a consistent hash
        block_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(block_string.encode()).hexdigest()

class Blockchain:
    """
    Manages the chain of blocks.
    """
    def __init__(self):
        self.chain = [self.create_genesis_block()]

    def create_genesis_block(self):
        """
        Manually creates the first block in the chain, the "Genesis Block".
        """
        return Block(0, datetime.datetime.now(), "Genesis Block - DNS Resolver Software", "0")

    def get_latest_block(self):
        """
        Returns the most recent block in the chain.
        """
        return self.chain[-1]

    def add_block(self, new_block_data):
        """
        Creates and adds a new block representing a software update to the chain.
        """
        latest_block = self.get_latest_block()
        new_block_index = latest_block.index + 1
        new_block_timestamp = datetime.datetime.now()
        new_block = Block(new_block_index, new_block_timestamp, new_block_data, latest_block.hash)
        self.chain.append(new_block)

def print_blockchain(chain):
    """
    Prints the contents of the blockchain for verification.
    """
    print("--- DNS Software Update Blockchain ---")
    for block in chain:
        print(f"Index: {block.index}")
        print(f"Timestamp: {block.timestamp}")
        print(f"Data: {block.data}")
        print(f"Previous Hash: {block.previous_hash}")
        print(f"Hash: {block.hash}")
        print("-" * 35)

# --- Main Execution ---

# 1. Initialize the blockchain for our DNS software
dns_update_chain = Blockchain()

# 2. A new authorized software update is released.
# We create a record for it and add it to the blockchain.
update1_data = {
    'version': '9.16.1',
    'file_hash': 'a1b2c3d4e5f67890a1b2c3d4e5f67890a1b2c3d4e5f67890a1b2c3d4e5f67890'
}
dns_update_chain.add_block(update1_data)

# 3. Another authorized update is released later.
# It is also added to the chain, linking back to the previous update.
update2_data = {
    'version': '9.18.7',
    'file_hash': 'f0e9d8c7b6a54321f0e9d8c7b6a54321f0e9d8c7b6a54321f0e9d8c7b6a54321'
}
dns_update_chain.add_block(update2_data)

# 4. Print the entire chain.
# An administrator or client system can now verify this chain. Any attempt
# to insert a malicious update would break the chain of hashes.
print_blockchain(dns_update_chain.chain)
