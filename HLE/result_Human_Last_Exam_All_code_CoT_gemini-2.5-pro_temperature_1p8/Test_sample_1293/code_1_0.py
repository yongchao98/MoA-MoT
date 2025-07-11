import hashlib
import json
import time

def simulate_blockchain_verification():
    """
    Simulates a secure software update mechanism for DNS software using a blockchain
    to prevent supply-chain attacks that would enable AitM.
    """
    # 1. The "Blockchain": A simulated immutable ledger. In reality, this would be a
    #    distributed and tamper-proof system. We use a list to represent the log of valid hashes.
    simulated_blockchain = []

    def create_software_package(version, content, is_tampered=False):
        """Helper function to create a software package dictionary."""
        if is_tampered:
            # An attacker injects malicious code to intercept/redirect DNS queries.
            content = content + "\n# MALICIOUS CODE: if query == 'bank.com' then redirect to 'evil.com'"
        package = {
            "component": "DNS-Resolver-Software",
            "version": version,
            "content": content,
            "timestamp": int(time.time()),
        }
        return package

    def calculate_secure_hash(package):
        """Calculates a SHA256 hash for a package to ensure its integrity."""
        # The package is converted to a consistent string format (sorted keys) before hashing.
        package_string = json.dumps(package, sort_keys=True).encode('utf-8')
        return hashlib.sha256(package_string).hexdigest()

    # 2. Developer Action: A legitimate developer releases a new version.
    print("--- Step 1: Legitimate Developer Action ---")
    legitimate_update = create_software_package(
        version="2.5.1",
        content="function handle_dns_query(query) { return securely_resolve(query); }"
    )
    legitimate_hash = calculate_secure_hash(legitimate_update)
    print(f"Developer created legitimate update version {legitimate_update['version']}.")
    print(f"The secure hash is: {legitimate_hash}")

    # The developer publishes this official hash to the public blockchain.
    simulated_blockchain.append(legitimate_hash)
    print("Developer has published this hash to the public blockchain.\n")

    # 3. Supply-Chain Attack: An adversary intercepts the update and creates a tampered version.
    #    This happens "in the middle" between the developer and the server.
    print("--- Step 2: Adversary Creates Tampered Update ---")
    tampered_update = create_software_package(
        version="2.5.1",
        content="function handle_dns_query(query) { return securely_resolve(query); }",
        is_tampered=True
    )
    tampered_hash = calculate_secure_hash(tampered_update)
    print(f"Adversary created a tampered package with malicious code.")
    print(f"Content: \"{tampered_update['content']}\"")
    print(f"The hash of the tampered package is: {tampered_hash}\n")


    # 4. Server Action: The DNS server's update mechanism verifies packages before applying them.
    def verify_and_apply_update(received_package, blockchain):
        """The server verifies a package against the blockchain before installation."""
        print(f"--- Verifying Package Version {received_package['version']} ---")
        hash_of_received_package = calculate_secure_hash(received_package)
        print(f"Calculated hash of received package: {hash_of_received_package}")

        # The crucial check: is the hash of the received software on the trusted blockchain?
        if hash_of_received_package in blockchain:
            print("  [SUCCESS] Verification passed: Hash found on the blockchain.")
            print("  Result: The update is authentic. Applying now.")
        else:
            print("  [DANGER] Verification FAILED: Hash NOT found on the blockchain.")
            print("  Result: Supply-chain attack detected! Rejecting the malicious update.")
        print("-" * 50)

    print("--- Step 3: Server's Automated Verification Process ---\n")
    # Scenario 1: Server receives the legitimate update.
    print("Scenario A: Verifying the legitimate package...")
    verify_and_apply_update(legitimate_update, simulated_blockchain)

    # Scenario 2: Server receives the tampered update from the attacker.
    print("\nScenario B: Verifying the tampered package...")
    verify_and_apply_update(tampered_update, simulated_blockchain)

if __name__ == "__main__":
    simulate_blockchain_verification()