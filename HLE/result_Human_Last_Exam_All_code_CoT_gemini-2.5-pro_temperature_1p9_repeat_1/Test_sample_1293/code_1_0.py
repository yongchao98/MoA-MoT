import hashlib

def verify_update_integrity():
    """
    Simulates the verification of a software update's integrity using
    cryptographic hashing, a core concept of the blockchain-based solution.
    """
    # This represents the official, authentic firmware data.
    legitimate_firmware_data = b'DNS_CLIENT_FIRMWARE_VERSION_3.1.4_BUILD_20240521'
    
    # The SHA-256 hash of the legitimate firmware is calculated and stored
    # on a secure, immutable ledger (e.g., a blockchain).
    # This hash is the "single source of truth".
    known_good_hash_from_ledger = hashlib.sha256(legitimate_firmware_data).hexdigest()
    
    print(f"[*] Secure Ledger contains the known-good hash: {known_good_hash_from_ledger}\n")

    # Scenario 1: The client receives the authentic update file.
    received_firmware_1 = legitimate_firmware_data
    print("--- SCENARIO 1: Verifying an authentic update ---")
    print(f"[*] Received firmware data: {received_firmware_1}")
    
    # The client calculates the hash of the received update.
    hash_of_received_firmware_1 = hashlib.sha256(received_firmware_1).hexdigest()
    print(f"[*] Calculated hash of received file: {hash_of_received_firmware_1}")

    # The client compares the calculated hash with the one from the secure ledger.
    if hash_of_received_firmware_1 == known_good_hash_from_ledger:
        print(">>> SUCCESS: Hash matches. The update is authentic and has not been tampered with.\n")
    else:
        print(">>> DANGER: HASH MISMATCH! The update is compromised. Aborting installation.\n")


    # Scenario 2: An adversary intercepts and modifies the update.
    tampered_firmware_data = b'DNS_CLIENT_FIRMWARE_VERSION_3.1.4_BUILD_20240521_MALICIOUS_IMPLANT'
    print("--- SCENARIO 2: Verifying a tampered update ---")
    print(f"[*] Received firmware data: {tampered_firmware_data}")

    # The client calculates the hash of the tampered update.
    hash_of_received_firmware_2 = hashlib.sha256(tampered_firmware_data).hexdigest()
    print(f"[*] Calculated hash of received file: {hash_of_received_firmware_2}")

    # The client compares the hash, and the mismatch reveals the tampering.
    if hash_of_received_firmware_2 == known_good_hash_from_ledger:
        print(">>> SUCCESS: Hash matches. The update is authentic and has not been tampered with.\n")
    else:
        print(">>> DANGER: HASH MISMATCH! The update is compromised. Aborting installation.\n")

if __name__ == '__main__':
    verify_update_integrity()