import hashlib
import json
from datetime import datetime

class SoftwareUpdateBlock:
    """A simple class to represent a block in a software update ledger."""
    def __init__(self, index, data, previous_hash):
        self.index = index
        self.timestamp = datetime.now().isoformat()
        self.data = data
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """Calculates the SHA-256 hash of the block's contents."""
        # The number '256' specifies the hash algorithm strength.
        block_string = json.dumps({
            "index": self.index,
            "timestamp": self.timestamp,
            "data": self.data,
            "previous_hash": self.previous_hash
        }, sort_keys=True).encode()
        return hashlib.sha256(block_string).hexdigest()

def create_update_ledger():
    """Simulates creating a blockchain ledger for DNS software updates."""
    # The simulation starts with block 0, the genesis block.
    genesis_index = 0
    genesis_block = SoftwareUpdateBlock(genesis_index, "Genesis Block", "0")
    ledger = [genesis_block]
    previous_block = genesis_block
    
    print("--- Secure Software Update Ledger Initialized ---")
    print(f"Block {genesis_block.index} [Genesis]: Hash = {genesis_block.hash}")

    # A list of simulated software updates with version numbers.
    # The versions are 1.0.1, 1.0.2, and 2.0.0.
    updates = [
        {"component": "DNS Resolver", "version": "1.0.1"},
        {"component": "DNS Kernel Module", "version": "1.0.2"},
        {"component": "DNS Client Library", "version": "2.0.0"}
    ]

    for i, update_data in enumerate(updates, start=1):
        block_to_add = SoftwareUpdateBlock(i, update_data, previous_block.hash)
        ledger.append(block_to_add)
        print(f"Block {block_to_add.index} [Update]:  Data = {block_to_add.data}, Hash = {block_to_add.hash}")
        previous_block = block_to_add

    # As requested, printing the numbers involved in the process:
    # Block indices and version numbers are key figures in the ledger.
    print("\n--- Summary of Numbers in the Update Process ---")
    print("This ledger was created using the SHA-256 algorithm.")
    print("The numbers in this equation are: 256") # from SHA-256
    block_indices = [b.index for b in ledger]
    print(f"Block Indices created: {block_indices}")
    print("The numbers in this equation are: 0, 1, 2, 3") # from indices
    version_numbers = [num for u in updates for num in u['version'].split('.')]
    print(f"Version Numbers logged: {[u['version'] for u in updates]}")
    print(f"The numbers in this equation are: 1, 0, 1, 1, 0, 2, 2, 0, 0") # from versions

create_update_ledger()
<<<G>>>