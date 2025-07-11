import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

print("Step-by-step reasoning for the diagnosis:")
print("1. Identify the trigger event: The patient's symptoms started after a 'difficult' colonoscopy.")
print("2. Analyze the key symptoms: The patient presented with a specific combination of Left Upper Quadrant (LUQ) pain and left-shoulder pain (Kehr's sign).")
print("3. Recognize signs of severe bleeding: There was a rapid development of tachycardia and hypotension (shock) and a significant drop in hemoglobin from 11.7 to 6.5 g/dL, indicating massive internal hemorrhage.")
print("4. Evaluate physical exam findings: Worsening LUQ tenderness and abdominal distension point to blood accumulating in the peritoneal cavity (hemoperitoneum).")
print("\nComparing the options:")
print("- A. Colonic perforation: Less likely. Doesn't explain the specific left shoulder pain or the purely hemorrhagic shock as well as other options.")
print("- B. Lower GI bleeding: Unlikely. This would present with bleeding from the rectum, not internal hemoperitoneum.")
print("- C. Splenic laceration: This is the most likely diagnosis. It is a known complication of difficult colonoscopies and perfectly explains the LUQ pain, Kehr's sign, and massive hemoperitoneum leading to shock.")
print("- D. Postpolypectomy syndrome: Ruled out because no polypectomy was performed.")
print("\nFinal conclusion: Splenic laceration is the diagnosis that fits all the clinical evidence.")
print("<<<C>>>")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output_str = output_buffer.getvalue()

# Print the captured output
print(output_str)