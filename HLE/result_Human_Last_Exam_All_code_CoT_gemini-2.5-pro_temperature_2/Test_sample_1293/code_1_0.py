import hashlib
import json
import time

def main():
    """
    Simulates using a blockchain to secure a software update process,
    preventing the installation of malicious code from a supply-chain attack.
    """

    # --- Part 1: Simulate the Blockchain and Software Publishing ---

    # A simple list represents our blockchain. In a real system, this would be
    # a distributed, cryptographically secured, and immutable ledger.
    blockchain = []

    def calculate_block_hash(index, timestamp, software_hash, previous_hash):
        """Calculates the hash of a block for chain integrity."""
        block_string = json.dumps({
            "index": index,
            "timestamp": timestamp,
            "software_hash": software_hash,
            "previous_hash": previous_hash
        }, sort_keys=True).encode()
        return hashlib.sha256(block_string).hexdigest()

    def hash_software(software_content):
        """Calculates the SHA-256 hash of a software package."""
        return hashlib.sha256(software_content.encode()).hexdigest()

    def publish_new_version(previous_block, software_content):
        """Hashes new software and adds its info as a new block to the chain."""
        index = previous_block["index"] + 1
        timestamp = time.time()
        sw_hash = hash_software(software_content)
        
        # In a real blockchain, the 'previous_hash' links blocks together securely.
        block_hash = calculate_block_hash(index, timestamp, sw_hash, previous_block["hash"])

        new_block = {
            "index": index,
            "timestamp": timestamp,
            "software_hash": sw_hash,  # The hash of the legitimate software update
            "previous_hash": previous_block["hash"],
            "hash": block_hash
        }
        
        blockchain.append(new_block)
        print(f"Published Block {index}:")
        print(f"  Legitimate Software Hash: {sw_hash}")
        print(f"  Block Hash for Ledger:    {new_block['hash']}\n")
        return new_block

    # --- Part 2: Simulate the Client-Side Verification ---

    def client_verify_and_install_update(software_content_received, trusted_blockchain):
        """A client verifies a received update against the trusted blockchain."""
        print("--- Client Verification Process ---")
        
        # 1. Get the latest official software hash from the trusted blockchain ledger.
        latest_block = trusted_blockchain[-1]
        official_hash = latest_block["software_hash"]
        print(f"Client retrieves latest official hash from blockchain: {official_hash}")

        # 2. Calculate the hash of the software update package it just received.
        received_hash = hash_software(software_content_received)
        print(f"Client calculates hash of the received software:       {received_hash}")
        
        # 3. Compare the hashes. If they match, the software is authentic.
        if received_hash == official_hash:
            print("\nResult: VERIFICATION SUCCESSFUL. Hashes match.")
            print("Action: Installing the authentic software update.")
        else:
            print("\nResult: VERIFICATION FAILED! Hashes DO NOT match.")
            print("Action: Rejecting the malicious or corrupted software update!")
        print("-" * 37 + "\n")

    # --- Main Simulation Execution ---

    # 1. Initialize the blockchain with a "Genesis Block" (the first block).
    genesis_block = {
        "index": 0, "timestamp": time.time(), "software_hash": "GENESIS_HASH", 
        "previous_hash": "0", "hash": calculate_block_hash(0, time.time(), "GENESIS_HASH", "0")
    }
    blockchain.append(genesis_block)
    print("Blockchain initialized.\n")

    # 2. Legitimate developer publishes a new DNS software version (v1.0).
    software_v1_content = "/* Legitimate DNS Server Code v1.0 */"
    last_block = publish_new_version(blockchain[-1], software_v1_content)

    # 3. Client A receives the legitimate update and verifies it.
    print("SCENARIO 1: Client receives the AUTHENTIC update.")
    client_verify_and_install_update(software_v1_content, blockchain)

    # 4. An attacker creates a malicious version and tries to push it to a client
    #    through a compromised download server (a supply-chain attack).
    malicious_software_content = "/* Malicious DNS Server Code v1.0 with implant */"

    # 5. Client B receives the malicious update and verifies it.
    print("SCENARIO 2: Client receives a MALICIOUS update.")
    client_verify_and_install_update(malicious_software_content, blockchain)

if __name__ == "__main__":
    main()