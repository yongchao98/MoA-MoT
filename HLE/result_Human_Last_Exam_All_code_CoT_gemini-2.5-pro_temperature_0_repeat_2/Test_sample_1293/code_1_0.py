import hashlib
import datetime

class Block:
    """
    A class to represent a block in our software update blockchain.
    """
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data  # e.g., {'version': '1.1.2', 'notes': 'Security fix...'}
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """
        Calculates the SHA-256 hash of the block's contents.
        """
        # We concatenate all block data into a single string to be hashed.
        block_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(block_string.encode()).hexdigest()

def create_genesis_block():
    """
    Creates the very first block in the chain.
    """
    return Block(0, datetime.datetime.now(), "Genesis Block: DNS Software v1.0.0 Initial Release", "0")

def create_new_block(previous_block, data):
    """
    Creates a new block to be added to the chain.
    """
    index = previous_block.index + 1
    timestamp = datetime.datetime.now()
    new_block = Block(index, timestamp, data, previous_block.hash)
    return new_block

# --- Simulation of Secure Update Process ---

# 1. Initialize the blockchain with the first-ever software release.
blockchain = [create_genesis_block()]
previous_block = blockchain[0]

# 2. A new security patch is released. It's added as a new block.
update_data_1 = {
    'software': 'DNS Client',
    'version': '1.0.1',
    'notes': 'Security patch for CVE-2023-12345.',
    'update_file_sha256': 'a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4e5f6a1b2'
}
block_to_add_1 = create_new_block(previous_block, update_data_1)
blockchain.append(block_to_add_1)
previous_block = block_to_add_1

# 3. A feature update is released. It's added as the next block.
update_data_2 = {
    'software': 'DNS Client',
    'version': '1.1.0',
    'notes': 'Added support for new record types.',
    'update_file_sha256': 'f6e5d4c3b2a1f6e5d4c3b2a1f6e5d4c3b2a1f6e5d4c3b2a1f6e5d4c3b2a1f6e5'
}
block_to_add_2 = create_new_block(previous_block, update_data_2)
blockchain.append(block_to_add_2)

# --- Output and Verification Demonstration ---
print("--- Secure Software Update Blockchain for DNS ---")
print("This simulation shows an immutable ledger for software updates.")
print("Each new update is a 'block' cryptographically chained to the previous one.\n")

for block in blockchain:
    print(f"Index: {block.index}")
    print(f"Timestamp: {block.timestamp}")
    print(f"Data: {block.data}")
    print(f"Previous Hash: {block.previous_hash}")
    print(f"Current Block Hash: {block.hash}\n")

# Fulfilling the requirement to show the "equation" for verification.
# A client verifies the latest block by recalculating its hash.
last_block = blockchain[-1]
print("--- Verification 'Equation' for the Final Block ---")
print("A client verifies integrity by checking if SHA256(block_contents) == stored_hash.\n")

# Printing each component of the final "equation"
index_num = last_block.index
timestamp_str = str(last_block.timestamp)
data_str = str(last_block.data)
prev_hash_str = last_block.previous_hash
final_hash_str = last_block.hash

# Re-calculate the hash for verification purposes
verification_string = str(index_num) + timestamp_str + data_str + prev_hash_str
calculated_hash = hashlib.sha256(verification_string.encode()).hexdigest()

print("Components for the final block's hash calculation:")
print(f"1. Index Number: {index_num}")
print(f"2. Timestamp String: '{timestamp_str}'")
print(f"3. Data String: '{data_str}'")
print(f"4. Previous Hash String: '{prev_hash_str}'")
print(f"\nResult of SHA256(concatenated components): '{calculated_hash}'")
print(f"Stored Hash in Block:                     '{final_hash_str}'")
print(f"\nVerification Successful: {calculated_hash == final_hash_str}")
