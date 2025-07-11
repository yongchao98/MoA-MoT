import hashlib

def get_software_hash(software_data):
    """Calculates the SHA-256 hash of software data."""
    return hashlib.sha256(software_data).hexdigest()

# --- Vendor Side ---
# The legitimate software update package.
original_software = b"This is the official DNS server software version 2.1"
# The vendor calculates the hash and publishes it to an immutable ledger (our simulated blockchain).
official_hash = get_software_hash(original_software)

# A simple dictionary to simulate the immutable blockchain ledger.
# In a real system, this would be a distributed and tamper-proof ledger.
blockchain_ledger = {
    "dns_server_v2.1": official_hash
}

print("--- Vendor publishes update hash to the ledger ---")
print(f"Software: dns_server_v2.1")
print(f"Official Hash: {blockchain_ledger['dns_server_v2.1']}\n")


# --- Attacker Action ---
# An adversary intercepts the update and injects a backdoor.
tampered_software = b"This is the official DNS server software version 2.1 ... with a malicious backdoor"


# --- Client Side Verification ---
print("--- Client attempts to verify downloaded software ---")

# Case 1: Client downloads the legitimate software
print("\nScenario 1: Verifying the legitimate software update...")
downloaded_software_1 = original_software
client_calculated_hash_1 = get_software_hash(downloaded_software_1)

# The "equation" is the comparison of the two hashes.
# We will print each "number" (hash) in the equation.
print(f"Hash calculated by client: {client_calculated_hash_1}")
print(f"Official hash from ledger: {blockchain_ledger['dns_server_v2.1']}")

if client_calculated_hash_1 == blockchain_ledger['dns_server_v2.1']:
    print("Result: VERIFIED. The hashes match. The software is authentic. It is safe to install.")
else:
    print("Result: REJECTED. The hashes DO NOT match. The software is compromised.")

# Case 2: Client downloads the tampered software
print("\nScenario 2: Verifying the tampered software update...")
downloaded_software_2 = tampered_software
client_calculated_hash_2 = get_software_hash(downloaded_software_2)

# Printing each "number" (hash) in the comparison equation.
print(f"Hash calculated by client: {client_calculated_hash_2}")
print(f"Official hash from ledger: {blockchain_ledger['dns_server_v2.1']}")

if client_calculated_hash_2 == blockchain_ledger['dns_server_v2.1']:
    print("Result: VERIFIED. The hashes match. The software is authentic. It is safe to install.")
else:
    print("Result: REJECTED. The hashes DO NOT match. The software is compromised and will not be installed.")
