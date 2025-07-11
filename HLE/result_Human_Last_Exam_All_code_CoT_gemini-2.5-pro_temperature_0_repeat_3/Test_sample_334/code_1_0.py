import hashlib
import os

class SecureBiometricModule:
    """
    A conceptual simulation of a FIPS 140-3 Level 5 compliant biometric module.
    This class demonstrates the multi-layered security checks required to address
    spoofing, coercion, and data exposure vulnerabilities.
    """
    def __init__(self, user_id, raw_biometric_data):
        print(f"--- Initializing Secure Module for user '{user_id}' ---")
        # In a real system, this would be a secure hardware enclave.
        self._secure_enclave = {}
        self._user_id = user_id
        self._is_revoked = False
        
        # Generate a salted hash of the biometric template.
        # This represents protecting data at rest.
        salt = os.urandom(16)
        self._secure_enclave['salt'] = salt
        self._secure_enclave['template_hash'] = self._hash_biometric(raw_biometric_data, salt)
        print("Biometric template securely hashed and stored in enclave.")
        print("-" * 20)

    def _hash_biometric(self, data, salt):
        """Simulates creating a secure hash of a biometric template."""
        # In reality, this would use fuzzy hashing and quantum-resistant algorithms.
        hasher = hashlib.sha256()
        hasher.update(salt)
        hasher.update(data.encode('utf-8'))
        return hasher.hexdigest()

    def revoke_credential(self):
        """Simulates revoking a compromised credential."""
        self._is_revoked = True
        # In a real FIPS module, this might trigger zeroization (wiping the key).
        self._secure_enclave = {}
        print(f"ALERT: Credentials for user '{self._user_id}' have been REVOKED.")

    def authenticate(self, biometric_input, liveness_check_passed, is_under_duress, context_signals):
        """
        Performs a full, multi-layered authentication check.
        
        Args:
            biometric_input (str): The live biometric data from the sensor.
            liveness_check_passed (bool): Result from liveness detection (e.g., 3D sensing).
            is_under_duress (bool): Result from coercion detection (e.g., abnormal heart rate, duress phrase).
            context_signals (dict): Contextual data like geolocation, device ID.
        """
        print(f"\n>>> Starting authentication attempt for '{self._user_id}'...")
        
        # 1. Check for Revocation (Addresses Data Exposure)
        if self._is_revoked:
            print("Step 1: Revocation Check... FAILED. Credential has been revoked.")
            print("Authentication FAILED.")
            return False
        print("Step 1: Revocation Check... PASSED.")

        # 2. Check for Coercion (Addresses Unwilling User)
        if is_under_duress:
            print("Step 2: Duress Check... FAILED. Signs of coercion detected.")
            # The system could silently alert security instead of just failing.
            print("Authentication FAILED. (Security silently alerted).")
            return False
        print("Step 2: Duress Check... PASSED.")

        # 3. Check for Liveness (Addresses Spoofing & Replay Attacks)
        if not liveness_check_passed:
            print("Step 3: Liveness Check... FAILED. Potential spoofing or replay attack detected.")
            print("Authentication FAILED.")
            return False
        print("Step 3: Liveness Check... PASSED.")

        # 4. Verify Biometric Data (Core Authentication)
        # This check only happens if all previous security layers have passed.
        print("Step 4: Verifying Biometric Template...")
        incoming_hash = self._hash_biometric(biometric_input, self._secure_enclave['salt'])
        if incoming_hash == self._secure_enclave['template_hash']:
            print("Biometric match... SUCCESS.")
            print("Authentication SUCCESSFUL.")
            return True
        else:
            print("Biometric match... FAILED.")
            print("Authentication FAILED.")
            return False

# --- Simulation Scenarios ---

# A unique physical trait for our user
user_fingerprint_data = "unique_fingerprint_pattern_12345"

# 1. Setup the secure module for a new user
module = SecureBiometricModule(user_id="alex", raw_biometric_data=user_fingerprint_data)

# 2. Scenario: Successful Authentication
print("\n--- SCENARIO 1: Normal, Successful Login ---")
module.authenticate(
    biometric_input="unique_fingerprint_pattern_12345", # Correct biometric
    liveness_check_passed=True, # Live person
    is_under_duress=False, # Not coerced
    context_signals={'location': 'home_wifi'}
)

# 3. Scenario: Spoofing Attack
print("\n--- SCENARIO 2: Spoofing Attack (e.g., 3D-printed finger) ---")
module.authenticate(
    biometric_input="unique_fingerprint_pattern_12345", # Data is correct
    liveness_check_passed=False, # But it's not from a live person
    is_under_duress=False,
    context_signals={'location': 'home_wifi'}
)

# 4. Scenario: Coercion Attack
print("\n--- SCENARIO 3: User Under Duress (Coercion) ---")
module.authenticate(
    biometric_input="unique_fingerprint_pattern_12345", # Correct biometric
    liveness_check_passed=True, # Live person
    is_under_duress=True, # But they are being forced
    context_signals={'location': 'unrecognized_van_moving_fast'}
)

# 5. Scenario: Stolen Data (Credential Revoked)
print("\n--- SCENARIO 4: Attempt with Stolen Data after Revocation ---")
# The database was breached, so the administrator revokes the credential.
module.revoke_credential()
# The attacker now tries to use the stolen data.
module.authenticate(
    biometric_input="unique_fingerprint_pattern_12345",
    liveness_check_passed=True,
    is_under_duress=False,
    context_signals={'location': 'attackers_laptop'}
)