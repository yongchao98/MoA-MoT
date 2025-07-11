import hashlib

def simulate_blockchain_verification():
    """
    Simulates verifying a DNS software update using a mock blockchain ledger.
    This demonstrates the core principle of ensuring software integrity to
    prevent supply-chain attacks that could enable AitM.
    """

    # A mock "blockchain ledger" storing the official, correct hashes for software versions.
    # In a real system, this would be a distributed, immutable ledger.
    blockchain_ledger = {
        "dns-server-v1.0.0": "d2f62b775402a799c8926987d6928e4e69b5523b105f7787265a88c4f6d3f2b1",
        "dns-server-v1.0.1": "1a4d8c7e0c4f8d2b9a7f6c5e3d2a1b0c9f8e7d6b5a4c3e2f1d0b9a8c7e6d5f4e"
    }

    # The legitimate software update binary (represented as bytes)
    software_version = "v1.0.1"
    legitimate_software_binary = b"This is the official DNS server software version 1.0.1 code."

    # A malicious binary, altered by an adversary in a supply-chain attack
    compromised_software_binary = b"This is the official DNS server software version 1.0.1 code. # MALICIOUS CODE INJECTED"

    print("--- Starting Software Integrity Verification ---\n")

    # --- Verification for the Legitimate Software ---
    print(f"Scenario 1: Verifying legitimate software update ({software_version})")

    # 1. Calculate the hash of the received software binary
    calculated_hash_legit = hashlib.sha256(legitimate_software_binary).hexdigest()

    # 2. Retrieve the official hash from the ledger
    official_hash = blockchain_ledger.get(f"dns-server-{software_version}", "Version not found in ledger")

    # This represents the "equation" for verification: Calculated_Hash == Official_Hash
    print(f"Calculated Hash: {calculated_hash_legit}")
    print(f"Official Hash (from Blockchain): {official_hash}")

    if calculated_hash_legit == official_hash:
        print("Result: VERIFICATION SUCCESSFUL. Software integrity is confirmed.\n")
    else:
        print("Result: VERIFICATION FAILED. Software has been TAMPERED with!\n")


    # --- Verification for the Compromised Software ---
    print(f"Scenario 2: Verifying compromised software update ({software_version})")

    # 1. Calculate the hash of the compromised binary
    calculated_hash_compromised = hashlib.sha256(compromised_software_binary).hexdigest()

    # 2. Retrieve the official hash (it's the same official hash)
    # The "equation" components for this scenario
    print(f"Calculated Hash: {calculated_hash_compromised}")
    print(f"Official Hash (from Blockchain): {official_hash}")

    if calculated_hash_compromised == official_hash:
        print("Result: VERIFICATION SUCCESSFUL. Software integrity is confirmed.\n")
    else:
        print("Result: VERIFICATION FAILED. Software has been TAMPERED with!\n")


if __name__ == '__main__':
    simulate_blockchain_verification()