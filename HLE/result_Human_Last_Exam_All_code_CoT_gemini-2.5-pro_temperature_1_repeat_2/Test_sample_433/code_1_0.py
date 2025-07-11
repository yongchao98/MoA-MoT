import datetime

def calculate_ampule_bud():
    """
    Explains and determines the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP <797> guidelines.
    """
    
    # Parameters from the user's question
    container_type = "Single-Dose Ampule"
    environment = "Sterile Environment (e.g., ISO Class 5)"
    
    # According to USP <797>, an ampule is an open-system container.
    # Once broken, it cannot be resealed to protect against contamination.
    # Therefore, it cannot be stored for any period.
    allowed_storage_hours = 0
    
    print("Analysis of Beyond-Use Date (BUD) for a Single-Dose Ampule")
    print("="*60)
    print(f"Container Type: {container_type}")
    print(f"Environment: {environment}")
    print("\nGuideline Summary (from USP <797>):")
    print("An opened ampule is an open system that cannot be resealed. To prevent microbial contamination, its contents must be used immediately. Any unused portion must be discarded and cannot be stored.")
    
    print("\nBUD Calculation:")
    # We represent "immediate use" as a storage time of 0 hours.
    puncture_time = datetime.datetime.now()
    bud_time = puncture_time + datetime.timedelta(hours=allowed_storage_hours)
    
    print(f"Allowed storage time after puncture: {allowed_storage_hours} hours")
    print(f"This means the contents must be used 'Immediately'.")
    
    print("\nFinal Conclusion:")
    print("The BUD for a single-dose container ampule from the time of puncture is Immediate Use.")
    print("Storage is not permitted.")


if __name__ == "__main__":
    calculate_ampule_bud()
