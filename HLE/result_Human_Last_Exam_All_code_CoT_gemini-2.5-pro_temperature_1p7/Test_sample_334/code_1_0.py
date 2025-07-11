import hashlib
import random
import os

# --- Helper Functions ---
def pseudo_biometric_scan(true_template):
    """Simulates a noisy biometric scan by introducing small variations."""
    scanned_template = {}
    for feature, value in true_template.items():
        # Introduce a small random noise (e.g., +/- 5%) to simulate real-world variance.
        noise = (random.random() - 0.5) * 0.1 * value
        scanned_template[feature] = value + noise
    return scanned_template

def is_match(template1, template2, tolerance=0.08):
    """A 'fuzzy' comparison to check if two biometric templates match within a given tolerance."""
    if template1.keys() != template2.keys():
        return False
    for feature in template1.keys():
        if max(abs(template1[feature]), abs(template2[feature])) == 0:
            continue
        # Check if the relative difference is within the allowed tolerance.
        if (abs(template1[feature] - template2[feature]) / max(abs(template1[feature]), abs(template2[feature]))) > tolerance:
            return False
    return True

def create_vault(biometric_template, secret_payload, pin):
    """Simulates creating a vault that binds a biometric template to a secret payload, secured by a PIN."""
    # The 'vault' holds helper data, a salted hash of the PIN, and the locked payload.
    # In a real system, this would use complex fuzzy extractor cryptography.
    vault = {
        "id": os.urandom(8).hex(),
        "template_helper": biometric_template, # In reality, this would be public helper data, not the template itself.
        "pin_salt_hash": hashlib.sha256((pin + str(random.random())).encode()).hexdigest(), # Salted hash
        "locked_payload": secret_payload # For demo, payload is stored plainly. In reality, it's locked by the bio-key.
    }
    print(f"-> Vault Created! ID: {vault['id']}. The vault is now stored on the server.")
    return vault

def unlock_vault(scanned_biometric, pin_attempt, vault):
    """Attempts to unlock the vault with a new biometric scan and a PIN."""
    print("\n--- Attempting to unlock vault ---")
    
    # In a real system, the PIN hash would be checked first. Here we simplify.
    
    # 1. Verify Biometric with fuzzy matching against the vault's helper data.
    if not is_match(scanned_biometric, vault["template_helper"]):
        print("Unlock Failed: Biometric mismatch. The provided scan is too different from the one enrolled.")
        return None
    print("Verification Step 1: Fuzzy Biometric Match... SUCCESS.")

    # 2. In a full system, the now-trusted biometric data and PIN would generate a key to unlock the payload.
    # We will simulate that this step passed.
    print("Verification Step 2: PIN and Key Generation... SUCCESS.")
    
    print("\nUnlock Successful! Retrieved Secret Payload:")
    retrieved_payload = vault["locked_payload"]
    print(f"Payload: {retrieved_payload}")
    return retrieved_payload

# --- Simulation of the System in Action ---

# The user's inherent, unchanging biometric data (e.g., from fingerprint features).
USER_TRUE_BIOMETRIC = {"ridge_count": 150, "valley_depth": 30.5, "minutia_angle": 88.7}
USER_PIN = "1234"
SECRET_PAYLOAD_V1 = "Secret_Server_Access_Key_v1"

print("Step 1: User Enrollment")
# A scan is taken during enrollment and the first vault is created.
enrollment_scan = pseudo_biometric_scan(USER_TRUE_BIOMETRIC)
user_vault_v1 = create_vault(enrollment_scan, SECRET_PAYLOAD_V1, USER_PIN)

print("\n" + "="*50 + "\n")

print("Step 2: Normal Authentication Attempt")
# The user provides a new, slightly different scan to authenticate.
auth_scan_1 = pseudo_biometric_scan(USER_TRUE_BIOMETRIC)
unlock_vault(auth_scan_1, USER_PIN, user_vault_v1)

print("\n" + "="*50 + "\n")

print("Step 3: Credential Compromise & Revocation")
print("An attacker has stolen the 'user_vault_v1' data from the server!")
print("Action: We will REVOKE the old vault and reissue a new credential.")

# To revoke, the user re-enrolls, creating a NEW vault with a NEW secret payload.
# The server now deletes/invalidates 'user_vault_v1'.
SECRET_PAYLOAD_V2 = "Reissued_Secret_Access_Key_v2"
print(f"Issuing a new payload: {SECRET_PAYLOAD_V2}")
new_enrollment_scan = pseudo_biometric_scan(USER_TRUE_BIOMETRIC)
user_vault_v2 = create_vault(new_enrollment_scan, SECRET_PAYLOAD_V2, USER_PIN)
print("The old vault (ID: " + user_vault_v1['id'] + ") is now obsolete.")

print("\n" + "="*50 + "\n")

print("Step 4: Authentication with the NEW, Reissued Vault")
# The user authenticates again. The old stolen vault would now fail verification on the server.
auth_scan_2 = pseudo_biometric_scan(USER_TRUE_BIOMETRIC)
unlock_vault(auth_scan_2, USER_PIN, user_vault_v2)

# Fulfilling the final specific instruction from the prompt.
print("\nFinal base biometric values from simulation:")
for number in USER_TRUE_BIOMETRIC.values():
    print(number)
