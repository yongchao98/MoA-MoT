import datetime

# --- Conditions ---
# The user is asking about the Beyond-Use Date (BUD) for a specific type of container.
# This is governed by pharmaceutical compounding standards (USP <797>).

container_type = "Single-Dose Container (Ampule)"
environment = "ISO Class 5 or cleaner (Sterile Environment)"

# --- Determination of BUD ---
# According to USP <797>, for a single-dose container opened in an ISO Class 5
# environment, the contents must be used within a specific time frame because
# it lacks antimicrobial preservatives.

bud_duration_hours = 12

# --- Output the Answer ---
print("Determining the Beyond-Use Date (BUD) based on USP <797> standards:")
print(f"Container Type: {container_type}")
print(f"Environment: {environment}")
print("-" * 50)
print("The standard states that after the initial puncture or opening, the container must be used within the following timeframe:")
print(f"BUD = {bud_duration_hours} hours")
print("-" * 50)
print(f"Therefore, the BUD for a single-dose ampule punctured in a sterile environment is 12 hours.")
