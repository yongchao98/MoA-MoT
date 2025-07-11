# Patient data with assigned priority scores based on clinical urgency.
# A lower score means a higher priority.
# Priority 1: Active neurologic deficit (surgical emergency).
# Priority 2: High spinal instability without neurologic deficit.
# Priority 3: Lesser spinal instability without neurologic deficit.

patients = [
    {
        'id': 1,
        'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.',
        'priority_score': 2
    },
    {
        'id': 2,
        'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
        'priority_score': 3
    },
    {
        'id': 3,
        'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.',
        'priority_score': 1
    }
]

# Sort patients by their priority score in ascending order (1 comes first).
sorted_patients = sorted(patients, key=lambda p: p['priority_score'])

# Print the final prioritization
print("Surgical Prioritization from top to lowest priority:")
priority_levels = ["Top priority", "Second priority", "Third priority"]

final_order = []
for i, patient in enumerate(sorted_patients):
    print(f"{priority_levels[i]}: Patient {patient['id']} ({patient['diagnosis']})")
    final_order.append(f"Patient {patient['id']}")

print("\nFinal Prioritization Order:")
print(', '.join(final_order))

# Determine the corresponding answer choice
answer_key = {
    'A': ['Patient 1', 'Patient 2', 'Patient 3'],
    'B': ['Patient 1', 'Patient 3', 'Patient 2'],
    'C': ['Patient 3', 'Patient 2', 'Patient 1'],
    'D': ['Patient 2', 'Patient 1', 'Patient 3'],
    'E': ['Patient 2', 'Patient 3', 'Patient 1'],
    'F': ['Patient 3', 'Patient 1', 'Patient 2']
}

correct_choice = ''
for choice, order in answer_key.items():
    if order == final_order:
        correct_choice = choice
        break

print(f"\nThis prioritization corresponds to answer choice {correct_choice}.")
print("The final answer is F. Patient 3 has a neurologic deficit and is the highest priority. Patient 1 has a highly unstable fracture and is the second priority. Patient 2 has a less unstable fracture and is the lowest priority.")
<<<F>>>