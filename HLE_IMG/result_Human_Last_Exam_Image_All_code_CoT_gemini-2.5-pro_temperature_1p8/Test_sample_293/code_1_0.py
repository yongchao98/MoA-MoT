def get_diagnosis():
    """
    Analyzes the key features from the histopathology image and provides a diagnosis.
    """
    diagnosis = "Crescentic glomerulonephritis"
    
    key_features = [
        "A large cellular crescent filling Bowman's space.",
        "Compression of the glomerular tuft.",
        "Hypercellularity within the glomerulus.",
        "Signs of fibrinoid necrosis (pink, amorphous material)."
    ]
    
    print(f"Histopathological Diagnosis: {diagnosis}")
    print("\nKey Observable Features:")
    for feature in key_features:
        print(f"- {feature}")

if __name__ == "__main__":
    get_diagnosis()