import hashlib
import os

def generate_file_hash(filename):
    """
    Generates a SHA-256 hash for a given file.
    This hash acts as a unique fingerprint to verify the file's integrity.
    In a blockchain system, this hash would be stored in an immutable ledger.
    """
    # Use a secure hash algorithm, SHA-256 is a standard choice.
    hasher = hashlib.sha256()
    
    # Open the file in binary read mode to ensure consistency across platforms
    try:
        with open(filename, 'rb') as f:
            # Read the file in chunks to handle large files efficiently
            while chunk := f.read(4096):
                hasher.update(chunk)
        return hasher.hexdigest()
    except FileNotFoundError:
        return None

def simulate_update_process():
    """
    Simulates the process of verifying a software update using a hash.
    """
    # --- STEP 1: On the Developer's Side (Create the update) ---
    software_filename = "dns_client_v2.zip"
    software_content = b"This is the binary content of the new DNS client software."
    
    # Create a dummy software file
    with open(software_filename, "wb") as f:
        f.write(software_content)
        
    # The developer calculates the hash of the official update file.
    # This hash is the "proof" of integrity.
    official_hash = generate_file_hash(software_filename)
    
    print("--- Software Update Verification Simulation ---")
    print(f"1. Developer creates software file: '{software_filename}'")
    print(f"2. Developer calculates the official SHA-256 hash.")
    print(f"   Official Hash: {official_hash}")
    print("   (This hash would be signed and stored on the blockchain)")
    print("-" * 50)

    # --- STEP 2: On the User's Side (Receive and verify the update) ---
    print("3. User's system downloads the software update file.")
    
    # The user's system calculates the hash of the downloaded file.
    downloaded_hash = generate_file_hash(software_filename)
    print(f"4. User's system calculates the hash of the downloaded file.")
    print(f"   Downloaded File Hash: {downloaded_hash}")
    
    # The user's system compares the calculated hash with the official one from the blockchain.
    print("5. System compares the local hash against the official hash from the blockchain.")
    
    if downloaded_hash == official_hash:
        print("\nResult: SUCCESS! The hash values match.")
        print("The software update is authentic and has not been tampered with.")
        print("The system can safely proceed with the installation.")
    else:
        print("\nResult: DANGER! The hash values DO NOT match.")
        print("The software update has been modified or corrupted.")
        print("The installation will be aborted to prevent a potential compromise.")
        
    # Clean up the dummy file
    os.remove(software_filename)

if __name__ == "__main__":
    simulate_update_process()
    print("\n<<<G>>>")
