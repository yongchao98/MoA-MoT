import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new stream
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output

print("The patient's presentation with sciatica-like pain in the L4-S1 distribution, which is exacerbated in certain positions but without remarkable X-ray findings, is highly suggestive of piriformis syndrome.")
print("Piriformis syndrome involves the compression or irritation of the sciatic nerve by the piriformis muscle.")
print("\nTo confirm this diagnosis through a provocative test, the physician aims to make the piriformis muscle contract to see if it reproduces the symptoms.")
print("The physical exam is being conducted with the patient in the left decubitus position (lying on the left side) with the right leg extended.")
print("\nThe primary action of the piriformis muscle when the hip is in an extended position is External Rotation.")
print("\nTherefore, the physician would ask the patient to externally rotate (turn outward) her extended right leg, and the physician would apply resistance to this movement.")
print("If this resisted action reproduces the patient's specific radiating pain, it helps to confirm the diagnosis of piriformis syndrome.")
print("\nLet's review the options:")
print("A. Abduction: Tests the gluteus medius/minimus.")
print("B. Adduction: Tests the adductor muscles.")
print("C. Internal Rotation: Opposite of the piriformis's primary action.")
print("D. External Rotation: The primary action of the piriformis in this position. Resisting this action will test the muscle.")
print("E. Flexion: Tests the iliopsoas and rectus femoris.")
print("F. Extension: Tests the gluteus maximus and hamstrings.")
print("\nThe correct action is External Rotation.")

# Get the content
output = captured_output.getvalue()
# Restore stdout
sys.stdout = original_stdout
# Print the captured output
print(output)
print("<<<D>>>")