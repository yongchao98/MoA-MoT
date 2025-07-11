import sys
import io

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("Step-by-step Deriviation of the Answer:")
print("1. Analyze Patient's Symptoms:")
print("   - Left-sided facial weakness (can't lift eyebrow): Points to Facial Nerve (CN VII) palsy.")
print("   - Loss of left acoustic reflex: Further supports CN VII palsy (affecting the stapedius muscle).")
print("   - Hoarseness & cough: Points to Vagus Nerve (CN X) palsy, specifically the recurrent laryngeal nerve which innervates the laryngeal muscles.")
print("2. Correlate with Imaging:")
print("   - A mass in the thoracic cavity can compress the left recurrent laryngeal nerve, which loops under the aorta in the chest. This provides a direct anatomical cause for the hoarseness.")
print("3. Evaluate the Anatomical Choices:")
print("   - A. Tensor tympani: Innervated by CN V. Less likely to be the primary cause of the acoustic reflex loss given the clear CN VII signs.")
print("   - B. Lateral rectus: Innervated by CN VI. Unrelated to the patient's symptoms.")
print("   - C. Intercostal muscles: Not related to the specific cranial nerve deficits.")
print("   - D. Cricothyroid: A key muscle of the larynx innervated by a branch of the Vagus Nerve (CN X). Dysfunction of the laryngeal muscles, including the cricothyroid, causes hoarseness. This structure directly links the patient's symptom (hoarseness) with the nerve (Vagus) being compressed by the thoracic mass.")
print("   - E. Stylopharyngeus: Innervated by CN IX. Unrelated to the primary symptoms.")
print("4. Conclusion:")
print("   - The cricothyroid muscle is the most relevant structure because it is directly involved in the symptom of hoarseness, which is explained by the compression of its nerve supply (the vagus nerve/recurrent laryngeal nerve) by the thoracic mass.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the user
print(output)
print("<<<D>>>")