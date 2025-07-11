import sys
from io import StringIO

# A class to redirect stdout
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

# Analysis of the clinical scenario
# The patient presents with complex psychiatric issues, polypharmacy, and an active substance use disorder (cannabis).
# The immediate priorities should be to address the most acute problems safely and gather necessary information.

# Step 1: Evaluate the options based on clinical priority.
# Option I: Counsel patient on stopping cannabis.
# Priority: High. The patient's heavy daily cannabis use is likely exacerbating his insomnia, anxiety, and blunting the effects of his medications. Addressing this is crucial.
option_I = "Counsel patient on stopping cannabis. (Priority: High)"

# Option II: Ask patient to request admission for detox off all psych meds.
# Priority: Low/Inappropriate. Abruptly stopping all medications is dangerous and can cause severe withdrawal and decompensation. Medication simplification is needed, but not a "detox".
option_II = "Ask patient to request admission for detox off all psych meds. (Priority: Low/Inappropriate)"

# Option III: Order a urine drug test.
# Priority: High. This is a standard of care to get objective data, confirm substance use, and screen for other unreported substances, especially given his history and request for stimulants.
option_III = "Order a urine drug test. (Priority: High)"

# Option IV: Prescribe melatonin for insomnia.
# Priority: High (Supportive). Insomnia is a major complaint. Melatonin is a safe, non-addictive option that can provide symptomatic relief and help build therapeutic rapport while addressing the root causes (like cannabis use).
option_IV = "Prescribe melatonin for insomnia. (Priority: High - Supportive)"

# Option V: Discontinue acamprosate and increase naltrexone.
# Priority: Low. There's no indication that his treatment for alcohol use disorder (which is in remission) needs to be changed. This is not an immediate priority.
option_V = "Discontinue acamprosate and increase naltrexone. (Priority: Low)"

# Option VI: Start atomoxetine.
# Priority: Invalid. The patient is already on atomoxetine 80 mg daily, as stated in the case description.
option_VI = "Start atomoxetine. (Priority: Invalid)"

# Step 2: Combine the three most appropriate, high-priority actions.
# The three most logical and immediate steps are I, III, and IV.
# I - Address the active substance use through counseling.
# III - Gather objective clinical data.
# IV - Provide safe, symptomatic relief for a chief complaint.

print("The three most appropriate and immediately prioritized treatment options are:")
print(f"1. Option {option_I}")
print(f"2. Option {option_III}")
print(f"3. Option {option_IV}")
print("\nThis combination addresses the patient's active substance use, gathers essential diagnostic information, and provides low-risk symptomatic relief for a major complaint.")
print("\nTherefore, the correct combination is I, III, and IV.")
print("\nFinal Answer Equation: I + III + IV")
print("<<<L>>>")