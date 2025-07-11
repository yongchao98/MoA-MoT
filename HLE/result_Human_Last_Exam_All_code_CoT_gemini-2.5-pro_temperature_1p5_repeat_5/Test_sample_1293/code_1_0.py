import hashlib
import datetime as date

class Block:
    """
    Represents a single block in our blockchain.
    Each block contains a software update ('data').
    """
    def __init__(self, index, timestamp, data, previous_hash):
        self.index = index
        self.timestamp = timestamp
        self.data = data
        self.previous_hash = previous_hash
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        """
        Calculates the SHA-256 hash of the block's contents.
        """
        hash_string = str(self.index) + str(self.timestamp) + str(self.data) + str(self.previous_hash)
        return hashlib.sha256(hash_string.encode()).hexdigest()

class Blockchain:
    """
    Manages the chain of blocks.
    """
    def __init__(self):
        # The chain starts with a "Genesis Block"
        self.chain = [self.create_genesis_block()]

    def create_genesis_block(self):
        """
        Manually creates the first block in the chain.
        """
        return Block(0, date.datetime.now(), "Genesis Block: DNS Software Initial Version 1.0", "0")

    def get_latest_block(self):
        """
        Returns the most recent block in the chain.
        """
        return self.chain[-1]

    def add_block(self, new_data):
        """
        Adds a new block (a new software update) to the chain.
        """
        latest_block = self.get_latest_block()
        new_block = Block(
            index=latest_block.index + 1,
            timestamp=date.datetime.now(),
            data=new_data,
            previous_hash=latest_block.hash
        )
        self.chain.append(new_block)
        print(f"Block #{new_block.index} has been added to the blockchain!")
        print(f"Data: {new_block.data}")
        print(f"Hash: {new_block.hash}\n")

    def is_chain_valid(self):
        """
        Verifies the integrity of the blockchain.
        It checks if each block's stored hash matches its calculated hash,
        and if it correctly points to the previous block's hash.
        """
        for i in range(1, len(self.chain)):
            current_block = self.chain[i]
            previous_block = self.chain[i-1]

            # Check if the block's own hash is valid
            if current_block.hash != current_block.calculate_hash():
                print(f"Chain is INVALID: Data in Block {current_block.index} has been tampered with.")
                return False
            
            # Check if the block points to the correct previous block hash
            if current_block.previous_hash != previous_block.hash:
                print(f"Chain is INVALID: The link between Block {previous_block.index} and Block {current_block.index} is broken.")
                return False
        
        print("Chain is VALID. All software updates are authentic.")
        return True


# --- Simulation ---

print("--- Creating a Secure DNS Software Update Chain ---\n")
dns_update_chain = Blockchain()

# Add legitimate software updates
print("Deploying new updates...\n")
dns_update_chain.add_block("Software Update: DNS Server v1.1 - Security Patch")
dns_update_chain.add_block("Software Update: DNS Server v1.2 - Performance Improvement")

print("--- Verifying Chain Integrity ---")
dns_update_chain.is_chain_valid()
print("\n" + "="*50 + "\n")


print("--- Simulating a Supply-Chain Attack ---")
print("An attacker attempts to tamper with a previous update...\n")

# Tampering with block 1
# This simulates an attacker changing a past update to inject malicious code.
# For example, changing a legitimate patch to a compromised one.
attacked_block = dns_update_chain.chain[1]
print(f"Original data in Block {attacked_block.index}: {attacked_block.data}")
attacked_block.data = "MALICIOUS UPDATE: DNS Server v1.1 - Trojan Implanted"
print(f"Tampered data in Block {attacked_block.index}: {attacked_block.data}\n")


print("--- Re-Verifying Chain Integrity After Attack ---")
# The validation will now fail because the hash of the tampered block no longer matches,
# breaking the cryptographic link to the next block.
dns_update_chain.is_chain_valid()
