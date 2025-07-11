import hashlib
import datetime

class Block:
    """
    Represents a single block in our blockchain.
    A block contains data for a software update.
    """
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """
        Calculates the hash of the block based on its contents.
        """
        # We create the string to be hashed, including all parts of the block
        # The f-string formatting implicitly converts non-string types to strings
        block_string = f"{self.index}{self.timestamp}{self.data}{self.previous_hash}"
        return hashlib.sha256(block_string.encode()).hexdigest()

def create_genesis_block():
    """
    Manually creates the first block in the chain.
    """
    return Block(0, datetime.datetime.now(), "Genesis Block: Initial DNS Software v1.0", "0")

def create_new_block(last_block, data):
    """
    Creates a new block to be added to the chain.
    """
    this_index = last_block.index + 1
    this_timestamp = datetime.datetime.now()
    this_hash = last_block.hash
    return Block(this_index, this_timestamp, data, this_hash)

def verify_chain(blockchain):
    """
    Verifies the integrity of the entire blockchain.
    """
    for i in range(1, len(blockchain)):
        current_block = blockchain[i]
        previous_block = blockchain[i-1]

        # 1. Verify that the stored hash of the current block is correct
        if current_block.hash != current_block.calculate_hash():
            print(f"\nVerification FAILED: Block {current_block.index} has been tampered with. Calculated hash does not match stored hash.")
            print(f"Stored hash:   {current_block.hash}")
            print(f"Recalculated:  {current_block.calculate_hash()}")
            return False

        # 2. Verify that the chain is linked correctly
        if current_block.previous_hash != previous_block.hash:
            print(f"\nVerification FAILED: Chain broken at Block {current_block.index}. The 'previous_hash' does not match the hash of Block {previous_block.index}.")
            print(f"Block {current_block.index} previous_hash: {current_block.previous_hash}")
            print(f"Block {previous_block.index} actual hash:   {previous_block.hash}")
            return False
            
    print("\nVerification PASSED: Blockchain integrity is intact.")
    return True

if __name__ == '__main__':
    # Create the blockchain and add the genesis block
    blockchain = [create_genesis_block()]
    previous_block = blockchain[0]

    print("--- Simulating Secure Software Updates using a Blockchain ---")
    print(f"Block #{previous_block.index} | Data: {previous_block.data} | Hash: {previous_block.hash}")

    # Add a new software update (patch)
    update_data_1 = "Software Update: DNS Security Patch v1.1"
    block_to_add_1 = create_new_block(previous_block, update_data_1)
    blockchain.append(block_to_add_1)
    print(f"Block #{block_to_add_1.index} | Data: {block_to_add_1.data} | Hash: {block_to_add_1.hash}")
    previous_block = block_to_add_1

    # Add another software update (new version)
    update_data_2 = "Software Update: DNS Major Release v2.0"
    block_to_add_2 = create_new_block(previous_block, update_data_2)
    blockchain.append(block_to_add_2)
    print(f"Block #{block_to_add_2.index} | Data: {block_to_add_2.data} | Hash: {block_to_add_2.hash}")
    
    # First, verify the pristine chain
    print("\n--- Verifying the original, untampered chain ---")
    verify_chain(blockchain)
    
    # Now, let's simulate a malicious actor tampering with an old update
    print("\n--- Tampering with Block 1 to simulate a supply-chain attack ---")
    blockchain[1].data = "MALICIOUS UPDATE: DNS Rootkit v1.1"
    print(f"Modified Block #1 Data: '{blockchain[1].data}'")
    
    # Verify the chain again to see if the tampering is detected
    print("\n--- Re-Verifying the tampered chain ---")
    verify_chain(blockchain)
<<<G>>>