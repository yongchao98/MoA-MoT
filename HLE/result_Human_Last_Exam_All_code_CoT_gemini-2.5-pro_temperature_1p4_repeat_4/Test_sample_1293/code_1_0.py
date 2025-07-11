import hashlib

class Block:
    """
    Represents a block in our secure update chain. Each block contains
    the software version, a hash of the update data, and the hash of the
    previous block, creating a cryptographic chain.
    """
    def __init__(self, version, data_hash, previous_hash):
        self.version = version
        self.data_hash = data_hash
        self.previous_hash = previous_hash
        self.block_hash = self.calculate_block_hash()

    def calculate_block_hash(self):
        """Calculates the hash for the current block."""
        # The block's integrity is based on its contents and its link to the past
        block_content = f"{self.version}{self.data_hash}{self.previous_hash}"
        return hashlib.sha256(block_content.encode()).hexdigest()

def create_update_chain():
    """
    Simulates creating a secure chain of software updates, printing each step
    to demonstrate the chain of trust.
    """
    # The first block, or "Genesis Block", for the initial software state
    # It has a "previous_hash" of "0" as it's the start of the chain.
    genesis_block = Block(
        version="1.0.0",
        data_hash=hashlib.sha256("Initial Trusted Software".encode()).hexdigest(),
        previous_hash="0"
    )
    
    update_chain = [genesis_block]
    print("--- Secure Software Update Ledger ---")
    print(f"Initializing Chain with Block 0 (Genesis Block)...")
    print(f"  Version: {genesis_block.version}")
    print(f"  Previous Hash: {genesis_block.previous_hash}")
    print(f"  Block Hash: {genesis_block.block_hash}\n")

    # A series of subsequent software updates
    updates = [
        {"version": "1.0.1", "data": "Security patch for DNS cache."},
        {"version": "1.1.0", "data": "Added support for new record types."},
        {"version": "1.1.1", "data": "Performance improvements and final integrity check."}
    ]

    for i, update in enumerate(updates):
        previous_block = update_chain[-1]
        
        # In a real system, this would be the hash of the actual update file
        current_data_hash = hashlib.sha256(update["data"].encode()).hexdigest()
        
        # Create the new block, linking it to the previous block's hash
        new_block = Block(
            version=update["version"],
            data_hash=current_data_hash,
            previous_hash=previous_block.block_hash
        )
        
        update_chain.append(new_block)
        print(f"Applying Update. Creating Block {i + 1}...")
        print(f"  Version: {new_block.version}")
        print(f"  Previous Hash (from Block {i}): {new_block.previous_hash}")
        print(f"  Block Hash: {new_block.block_hash}\n")
        
    # Explain the final verification process
    final_block = update_chain[-1]
    print("--- Final System Integrity Verification ---")
    print(f"To trust the final version ({final_block.version}), we verify its hash.")
    print("This hash is derived from the update's content AND the hash of the entire history before it.")
    
    # Show the final "equation" for the hash calculation
    print("\nFinal Hash Calculation:")
    print(f"HASH = SHA256( Version + DataHash + PreviousBlockHash )")
    print(f"HASH = SHA256( '{final_block.version}' + '{final_block.data_hash}' + '{final_block.previous_hash}' )")
    print(f"Final Calculated Hash = {final_block.block_hash}")
    print("\nAn implant from a tampered update would fail this verification, as the hashes would not match.")

if __name__ == '__main__':
    create_update_chain()
